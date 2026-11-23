/* Radio visibility engine for the browser (runs in a Web Worker; also loadable in Node for tests).
 * JavaScript port of radio_los.py: same models, same defaults.
 *  - line of sight + 1st Fresnel zone clearance over terrain + buildings, 4/3 earth
 *  - ITU-R P.526 knife-edge diffraction, Deygout method (3 levels)
 *  - Okumura-Hata urban model for high site <-> street-level links
 *  - repeater site search over rooftops / high ground
 */
'use strict';

const R_EFF = 6371000 * 4 / 3;
const STEP = 5.0;

const DEFAULTS = { freq: 446, tx_dbm: 37, sens_dbm: -118, gain_tx: 0, gain_rx: 0, losses: 3,
  repeaters: false, rep_tx_dbm: 44, rep_gain: 3, rep_mast: 3, min_roof: 25, serve_margin: 10,
  search_pad: 2000, max_candidates: 600 };
const LIMITS = { freq: [30, 6000], tx_dbm: [-10, 60], sens_dbm: [-150, -60], gain_tx: [-10, 30], gain_rx: [-10, 30],
  losses: [0, 40], rep_tx_dbm: [0, 60], rep_gain: [-10, 30], rep_mast: [0, 100], min_roof: [0, 300],
  serve_margin: [-20, 60], search_pad: [0, 10000], max_candidates: [10, 3000] };

/* ---------- UTM (WGS84) via Krüger series, accurate to well below 1 mm in the zone ---------- */
const UTM = (() => {
  const a = 6378137, f = 1 / 298.257223563, k0 = 0.9996;
  const n = f / (2 - f), n2 = n * n, n3 = n2 * n, n4 = n3 * n, n5 = n4 * n, n6 = n5 * n;
  const A = a / (1 + n) * (1 + n2 / 4 + n4 / 64 + n6 / 256);
  const al = [0, n / 2 - 2 * n2 / 3 + 5 * n3 / 16 + 41 * n4 / 180 - 127 * n5 / 288 + 7891 * n6 / 37800,
    13 * n2 / 48 - 3 * n3 / 5 + 557 * n4 / 1440 + 281 * n5 / 630 - 1983433 * n6 / 1935360,
    61 * n3 / 240 - 103 * n4 / 140 + 15061 * n5 / 26880 + 167603 * n6 / 181440,
    49561 * n4 / 161280 - 179 * n5 / 168 + 6601661 * n6 / 7257600,
    34729 * n5 / 80640 - 3418889 * n6 / 1995840, 212378941 * n6 / 319334400];
  const be = [0, n / 2 - 2 * n2 / 3 + 37 * n3 / 96 - n4 / 360 - 81 * n5 / 512 + 96199 * n6 / 604800,
    n2 / 48 + n3 / 15 - 437 * n4 / 1440 + 46 * n5 / 105 - 1118711 * n6 / 3870720,
    17 * n3 / 480 - 37 * n4 / 840 - 209 * n5 / 4480 + 5569 * n6 / 90720,
    4397 * n4 / 161280 - 11 * n5 / 504 - 830251 * n6 / 7257600,
    4583 * n5 / 161280 - 108847 * n6 / 3991680, 20648693 * n6 / 638668800];
  const e = Math.sqrt(f * (2 - f));
  return {
    fwd(lon, lat, zone) {
      const lon0 = (zone * 6 - 183) * Math.PI / 180;
      const phi = lat * Math.PI / 180, lam = lon * Math.PI / 180 - lon0;
      const t = Math.sinh(Math.atanh(Math.sin(phi)) - e * Math.atanh(e * Math.sin(phi)));
      const xi = Math.atan2(t, Math.cos(lam)), eta = Math.atanh(Math.sin(lam) / Math.sqrt(1 + t * t));
      let X = xi, Y = eta;
      for (let j = 1; j <= 6; j++) { X += al[j] * Math.sin(2 * j * xi) * Math.cosh(2 * j * eta); Y += al[j] * Math.cos(2 * j * xi) * Math.sinh(2 * j * eta); }
      return [500000 + k0 * A * Y, k0 * A * X];
    },
    inv(x, y, zone) {
      const lon0 = (zone * 6 - 183) * Math.PI / 180;
      const xi = y / (k0 * A), eta = (x - 500000) / (k0 * A);
      let X = xi, Y = eta;
      for (let j = 1; j <= 6; j++) { X -= be[j] * Math.sin(2 * j * xi) * Math.cosh(2 * j * eta); Y -= be[j] * Math.cos(2 * j * xi) * Math.sinh(2 * j * eta); }
      const tau0 = Math.sin(X) / Math.sqrt(Math.sinh(Y) ** 2 + Math.cos(X) ** 2);
      let tau = tau0;
      for (let i = 0; i < 6; i++) {   // Newton iteration for conformal -> geodetic latitude
        const s = Math.sinh(e * Math.atanh(e * tau / Math.sqrt(1 + tau * tau)));
        const ti = tau * Math.sqrt(1 + s * s) - s * Math.sqrt(1 + tau * tau);
        const d = (tau0 - ti) / Math.sqrt(1 + ti * ti) * (1 + (1 - e * e) * tau * tau) / ((1 - e * e) * Math.sqrt(1 + tau * tau));
        tau += d; if (Math.abs(d) < 1e-14) break;
      }
      return [(lon0 + Math.atan2(Math.sinh(Y), Math.cos(X))) * 180 / Math.PI, Math.atan(tau) * 180 / Math.PI];
    }
  };
})();

/* ---------- surface model ---------- */
class Model {
  constructor(meta, surfaceDm, groundDm, taggedBits) {
    Object.assign(this, { meta, S: surfaceDm, G: groundDm, T: taggedBits });
    this.W = meta.width; this.H = meta.height; this.x0 = meta.x0; this.y0 = meta.y0; this.res = meta.res;
    this.gf = meta.ground_factor; this.GW = meta.ground_width; this.GH = meta.ground_height; this.zone = meta.utm_zone;
    const c = [[0, 0], [this.W, 0], [0, this.H], [this.W, this.H]].map(([cc, rr]) => this.toLL(this.x0 + cc * this.res, this.y0 - rr * this.res));
    this.bounds = [[Math.max(c[2][1], c[3][1]), Math.max(c[0][0], c[2][0])], [Math.min(c[0][1], c[1][1]), Math.min(c[1][0], c[3][0])]];
  }
  toXY(lon, lat) { return UTM.fwd(lon, lat, this.zone); }
  toLL(x, y) { return UTM.inv(x, y, this.zone); }
  col(x) { return (x - this.x0) / this.res; }
  row(y) { return (this.y0 - y) / this.res; }
  inside(x, y) { const r = this.row(y), c = this.col(x); return r >= 0 && r < this.H && c >= 0 && c < this.W; }
  idx(x, y) {
    let r = Math.round(this.row(y)), c = Math.round(this.col(x));
    r = r < 0 ? 0 : r >= this.H ? this.H - 1 : r; c = c < 0 ? 0 : c >= this.W ? this.W - 1 : c;
    return r * this.W + c;
  }
  surf(x, y) { return this.S[this.idx(x, y)] * 0.1; }
  tagged(x, y) { const i = this.idx(x, y); return (this.T[i >> 3] >> (i & 7)) & 1; }
  ground(x, y) {   // bilinear on the 30 m grid (cell centres)
    let u = this.col(x) / this.gf - 0.5, v = this.row(y) / this.gf - 0.5;
    u = Math.min(Math.max(u, 0), this.GW - 1.001); v = Math.min(Math.max(v, 0), this.GH - 1.001);
    const i = Math.floor(u), j = Math.floor(v), fu = u - i, fv = v - j, G = this.G, W = this.GW;
    const a = G[j * W + i], b = G[j * W + i + 1], c = G[(j + 1) * W + i], d = G[(j + 1) * W + i + 1];
    return 0.1 * ((a * (1 - fu) + b * fu) * (1 - fv) + (c * (1 - fu) + d * fu) * fv);
  }
}

/* ---------- propagation ---------- */
const J = v => v > -0.78 ? 6.9 + 20 * Math.log10(Math.sqrt((v - 0.1) ** 2 + 1) + v - 0.1) : 0;

function deygout(d, z, i0, i1, ha, hb, lam, depth) {
  if (i1 - i0 + 1 < 3 || depth > 2) return 0;
  const D = d[i1] - d[i0];
  let best = -Infinity, bi = -1;
  for (let i = i0 + 1; i < i1; i++) {
    const d1 = d[i] - d[i0], d2 = D - d1;
    const v = (z[i] - (ha + (hb - ha) * d1 / D)) * Math.sqrt(2 * D / (lam * d1 * d2));
    if (v > best) { best = v; bi = i; }
  }
  if (best <= -0.78) return 0;
  return J(best) + deygout(d, z, i0, bi, ha, z[bi], lam, depth + 1) + deygout(d, z, bi, i1, z[bi], hb, lam, depth + 1);
}

function hataUrban(f, hb, hm, dkm) {
  hb = Math.min(Math.max(hb, 30), 200); hm = Math.min(Math.max(hm, 1), 10); dkm = Math.max(dkm, 1);
  const lf = Math.log10(f), ahm = (1.1 * lf - 0.7) * hm - (1.56 * lf - 0.8);
  return 69.55 + 26.16 * lf - 13.82 * Math.log10(hb) - ahm + (44.9 - 6.55 * Math.log10(hb)) * Math.log10(dkm);
}

function link(M, a, b, freq, exR = 10) {
  const lam = 299.792458 / freq;
  const D = Math.hypot(b.x - a.x, b.y - a.y);
  const n = Math.max(Math.floor(D / STEP), 2), N = n + 1;
  const d = new Float64Array(N), z = new Float64Array(N);
  let worst = Infinity, wi = -1;
  for (let i = 0; i < N; i++) {
    const t = i / n, x = a.x + (b.x - a.x) * t, y = a.y + (b.y - a.y) * t, di = t * D;
    d[i] = di;
    const mid = di > exR && di < D - exR;
    z[i] = (mid ? M.surf(x, y) : -1e4) + di * (D - di) / (2 * R_EFF);
    if (mid) {
      const f1 = Math.sqrt(lam * di * (D - di) / Math.max(D, 1));
      const clr = f1 > 0 ? (a.z + (b.z - a.z) * t - z[i]) / f1 : Infinity;
      if (clr < worst) { worst = clr; wi = i; }
    }
  }
  let tagged = true, building = false;
  if (wi >= 0) {
    const t = wi / n, x = a.x + (b.x - a.x) * t, y = a.y + (b.y - a.y) * t;
    tagged = !!M.tagged(x, y); building = M.surf(x, y) - M.ground(x, y) > 2;
  } else worst = 99;
  const fspl = 20 * Math.log10(Math.max(D, 1) / 1000) + 20 * Math.log10(freq) + 32.44;
  const diff = deygout(d, z, 0, N - 1, a.z, b.z, lam, 0);
  let loss = fspl + diff, model = 'terrain+buildings';
  const ga = a.z - M.ground(a.x, a.y), gb = b.z - M.ground(b.x, b.y);
  const hi = Math.max(ga, gb), lo = Math.min(ga, gb);
  if (hi >= 30 && lo <= 10 && diff > 0) {
    const hl = hataUrban(freq, hi, lo, D / 1000);
    if (hl < loss) { loss = Math.max(hl, fspl); model = 'urban (Hata)'; }
  }
  return { dist_km: D / 1000, los: worst > 0, fresnel: worst, fspl, diff, loss, model, tagged, building };
}

const classify = m => m >= 20 ? 'good' : m >= 10 ? 'marginal' : m >= 0 ? 'weak' : 'no link';
const r1 = v => Math.round(v * 10) / 10;

function siteType(r) {
  if (!r.roof) return 'ground mast';
  return r.bh >= 190 ? 'check: probably a tower/mast, not a roof' : 'rooftop';
}

/* ---------- API ---------- */
function params(p) {
  const out = Object.assign({}, DEFAULTS);
  for (const [k, v0] of Object.entries(p || {})) {
    if (k === 'repeaters') { out[k] = !!v0; continue; }
    if (!(k in LIMITS) || v0 === null || v0 === '') continue;
    const v = Number(v0), [lo, hi] = LIMITS[k];
    if (!Number.isFinite(v) || v < lo || v > hi) throw new Error(`${k} = ${v0} is outside ${lo}..${hi}`);
    out[k] = k === 'max_candidates' ? Math.round(v) : v;
  }
  return out;
}

function placePoints(M, raw) {
  const names = new Set();
  return raw.map((p, i) => {
    let name = String(p.name || `P${i + 1}`).trim().slice(0, 40);
    if (names.has(name)) name = `${name} (${i + 1})`;
    names.add(name);
    const q = { name, lat: +p.lat, lon: +p.lon, h: p.h === undefined ? 1.5 : +p.h, roof: !!p.roof };
    if (!Number.isFinite(q.lat) || !Number.isFinite(q.lon)) throw new Error(`${name}: bad coordinates`);
    if (!(q.h >= 0 && q.h <= 500)) throw new Error(`${name}: antenna height ${p.h} m is not plausible`);
    const [x, y] = M.toXY(q.lon, q.lat);
    if (!M.inside(x, y)) throw new Error(`${name} is outside the map area of the surface model`);
    const base = q.roof ? M.surf(x, y) : M.ground(x, y);
    return Object.assign(q, { x, y, base, z: base + q.h });
  });
}

function computeLinks(M, pts, P) {
  const budget = P.tx_dbm + P.gain_tx + P.gain_rx - P.losses, rows = [];
  for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++) {
    const L = link(M, pts[i], pts[j], P.freq), rx = budget - L.loss, margin = rx - P.sens_dbm;
    rows.push({ a: pts[i].name, b: pts[j].name, dist_km: Math.round(L.dist_km * 100) / 100, los: L.los,
      fresnel: Math.round(L.fresnel * 100) / 100, fspl: r1(L.fspl), diff: r1(L.diff), model: L.model,
      rx_dbm: r1(rx), margin: r1(margin), status: classify(margin),
      obstacle: L.building ? 'building' : 'terrain', obstacle_height: L.tagged ? 'tagged' : 'estimated' });
  }
  return rows;
}

function repeaterSearch(M, pts, P, onProgress, top = 15) {
  const xs = pts.map(p => p.x), ys = pts.map(p => p.y), pad = P.search_pad;
  const r0 = Math.max(Math.floor(M.row(Math.max(...ys) + pad)), 0), c0 = Math.max(Math.floor(M.col(Math.min(...xs) - pad)), 0);
  const r1_ = Math.min(Math.floor(M.row(Math.min(...ys) - pad)), M.H - 1), c1 = Math.min(Math.floor(M.col(Math.max(...xs) + pad)), M.W - 1);
  const h = r1_ - r0, w = c1 - c0;
  if (h < 3 || w < 3) return [];
  const S = new Float32Array(h * w), hb = new Float32Array(h * w);
  for (let r = 0; r < h; r++) for (let c = 0; c < w; c++) {
    const x = M.x0 + (c0 + c) * M.res, y = M.y0 - (r0 + r) * M.res, s = M.S[(r0 + r) * M.W + c0 + c] * 0.1;
    S[r * w + c] = s; hb[r * w + c] = s - M.ground(x, y);
  }
  // 15x15 maximum filter (separable, edges clamped)
  const K = 7, tmp = new Float32Array(h * w), mx = new Float32Array(h * w);
  for (let r = 0; r < h; r++) for (let c = 0; c < w; c++) {
    let m = -Infinity; for (let k = Math.max(c - K, 0); k <= Math.min(c + K, w - 1); k++) m = Math.max(m, S[r * w + k]); tmp[r * w + c] = m;
  }
  for (let c = 0; c < w; c++) for (let r = 0; r < h; r++) {
    let m = -Infinity; for (let k = Math.max(r - K, 0); k <= Math.min(r + K, h - 1); k++) m = Math.max(m, tmp[k * w + c]); mx[r * w + c] = m;
  }
  const cands = [];
  for (let r = 0; r < h; r++) for (let c = 0; c < w; c++) {
    const i = r * w + c;
    if (S[i] !== mx[i]) continue;
    const gi = (r0 + r) * M.W + c0 + c, tg = (M.T[gi >> 3] >> (gi & 7)) & 1;
    if ((hb[i] >= P.min_roof && tg) || hb[i] < 2) cands.push([S[i], r, c]);
  }
  cands.sort((p, q) => q[0] - p[0]);
  const up = P.tx_dbm + P.gain_tx + P.rep_gain - P.losses, dn = P.rep_tx_dbm + P.rep_gain + P.gain_rx - P.losses;
  const budget = Math.min(up, dn), res = [], list = cands.slice(0, P.max_candidates);
  list.forEach(([s, r, c], k) => {
    const x = M.x0 + (c0 + c + 0.5) * M.res, y = M.y0 - (r0 + r + 0.5) * M.res;
    const surfC = M.S[(r0 + r) * M.W + c0 + c] * 0.1, bh = hb[r * w + c];
    const cd = { x, y, z: surfC + P.rep_mast, base: surfC, roof: bh >= 2, bh };
    const margins = pts.map(p => budget - link(M, cd, p, P.freq).loss - P.sens_dbm);
    res.push(Object.assign(cd, { margins, served: margins.filter(m => m >= P.serve_margin).length, worst: Math.min(...margins) }));
    if (onProgress && k % 25 === 0) onProgress(k / list.length);
  });
  res.sort((p, q) => q.served - p.served || q.worst - p.worst);
  const picked = [];
  for (const r of res) {
    if (picked.every(q => Math.hypot(r.x - q.x, r.y - q.y) > 300)) picked.push(r);
    if (picked.length >= top) break;
  }
  return picked.map((r, i) => {
    const [lon, lat] = M.toLL(r.x, r.y);
    return { rank: i + 1, lat: Math.round(lat * 1e6) / 1e6, lon: Math.round(lon * 1e6) / 1e6, type: siteType(r),
      building_m: Math.round(r.bh), site_asl: Math.round(r.base), served: r.served, total: pts.length, worst: r1(r.worst),
      margins: Object.fromEntries(pts.map((p, j) => [p.name, r1(r.margins[j])])) };
  });
}

function analyze(M, body, onProgress) {
  const P = params(body.params), pts = placePoints(M, body.points || []);
  if (pts.length < 1) throw new Error('Add at least one point');
  const links = pts.length >= 2 ? computeLinks(M, pts, P) : [];
  const repeaters = P.repeaters ? repeaterSearch(M, pts, P, onProgress) : [];
  return { points: pts.map(p => ({ name: p.name, lat: p.lat, lon: p.lon, h: p.h, roof: p.roof, base: r1(p.base), antenna_asl: r1(p.z) })),
    links, repeaters, params: P };
}

function profile(M, body) {
  const P = params(body.params), [a, b] = placePoints(M, [body.a, body.b]), n = 400;
  const lam = 299.792458 / P.freq, D = Math.hypot(b.x - a.x, b.y - a.y);
  const seg = Math.max(Math.floor(D / n / STEP), 1), half = Math.floor(seg / 2);
  const out = { d: [], ground: [], surface: [], line: [], f_lo: [], f_hi: [] };
  for (let k = 0; k <= n; k++) {
    const t = k / n, dk = t * D, bulge = dk * (D - dk) / (2 * R_EFF);
    let s = -Infinity;
    for (let m = Math.max(k * seg - half, 0); m <= Math.min(k * seg + half, n * seg); m++) {
      const tt = m / (n * seg); s = Math.max(s, M.surf(a.x + (b.x - a.x) * tt, a.y + (b.y - a.y) * tt));
    }
    const line = a.z + (b.z - a.z) * t, f1 = Math.sqrt(lam * dk * (D - dk) / Math.max(D, 1));
    out.d.push(r1(dk)); out.ground.push(r1(M.ground(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t) + bulge));
    out.surface.push(r1(s + bulge)); out.line.push(r1(line)); out.f_lo.push(r1(line - f1)); out.f_hi.push(r1(line + f1));
  }
  const L = link(M, a, b, P.freq);
  return Object.assign(out, { za: r1(a.z), zb: r1(b.z), names: [a.name, b.name], los: L.los,
    fresnel: Math.round(L.fresnel * 100) / 100, dist_km: Math.round(L.dist_km * 100) / 100, freq: P.freq });
}

function elevation(M, lat, lon) {
  const [x, y] = M.toXY(lon, lat);
  if (!M.inside(x, y)) return { inside: false };
  const g = M.ground(x, y), s = M.surf(x, y);
  return { inside: true, ground: r1(g), surface: r1(s), building: s - g > 1 ? r1(s - g) : 0, height_tagged: !!M.tagged(x, y) };
}

/* ---------- data loading ---------- */
async function gunzipMaybe(buf) {
  const u = new Uint8Array(buf);
  if (u[0] !== 0x1f || u[1] !== 0x8b) return buf;   // server already decoded it
  const stream = new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'));
  return await new Response(stream).arrayBuffer();
}
function undelta(buf, w, h) {
  const a = new Int16Array(buf);
  if (a.length !== w * h) throw new Error('model file has unexpected size');
  for (let r = 0; r < h; r++) { let acc = 0; const o = r * w; for (let c = 0; c < w; c++) { acc += a[o + c]; a[o + c] = acc; } }
  return a;
}
async function fetchBuf(url, onBytes) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Cannot download ${url} (${res.status})`);
  if (!res.body || !onBytes) return await res.arrayBuffer();
  const reader = res.body.getReader(), chunks = []; let got = 0;
  for (;;) { const { done, value } = await reader.read(); if (done) break; chunks.push(value); got += value.length; onBytes(value.length); }
  const out = new Uint8Array(got); let o = 0; for (const c of chunks) { out.set(c, o); o += c.length; }
  return out.buffer;
}
async function loadModel(base, onProgress) {
  const meta = await (await fetch(base + 'meta.json')).json();
  const total = Object.values(meta.bytes).reduce((s, v) => s + v, 0); let got = 0;
  const tick = n => { got += n; onProgress && onProgress(Math.min(got / total, 1)); };
  const [s, g, t] = await Promise.all(['surface', 'ground', 'tagged'].map(k => fetchBuf(base + meta.files[k], tick)));
  const S = undelta(await gunzipMaybe(s), meta.width, meta.height);
  const G = undelta(await gunzipMaybe(g), meta.ground_width, meta.ground_height);
  const T = new Uint8Array(await gunzipMaybe(t));
  return new Model(meta, S, G, T);
}

/* ---------- worker glue ---------- */
if (typeof self !== 'undefined' && typeof importScripts === 'function') {
  let M = null;
  self.onmessage = async ev => {
    const { id, type, body } = ev.data;
    const reply = (ok, data) => self.postMessage({ id, ok, data });
    try {
      if (type === 'init') {
        M = await loadModel(body.base, p => self.postMessage({ type: 'progress', what: 'download', p }));
        reply(true, { bounds: M.bounds, defaults: DEFAULTS });
      } else if (!M) throw new Error('Model not loaded yet');
      else if (type === 'analyze') reply(true, analyze(M, body, p => self.postMessage({ type: 'progress', what: 'repeaters', p })));
      else if (type === 'profile') reply(true, profile(M, body));
      else if (type === 'elevation') reply(true, elevation(M, +body.lat, +body.lon));
      else throw new Error('unknown request ' + type);
    } catch (e) { reply(false, e.message || String(e)); }
  };
}
if (typeof module !== 'undefined') module.exports = { UTM, Model, analyze, profile, elevation, loadModel, gunzipMaybe, undelta, DEFAULTS };
