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
  search_pad: 2000, max_candidates: 600, wall_extra: 3, wall_er: 5, l_urban: 6.8, q90: 0.5, q_nu: 1.5, x_br: 300 };
const LIMITS = { freq: [30, 6000], tx_dbm: [-10, 60], sens_dbm: [-150, -60], gain_tx: [-10, 30], gain_rx: [-10, 30],
  losses: [0, 40], rep_tx_dbm: [0, 60], rep_gain: [-10, 30], rep_mast: [0, 100], min_roof: [0, 300],
  serve_margin: [-20, 60], search_pad: [0, 10000], max_candidates: [10, 3000], wall_extra: [0, 30], wall_er: [1.5, 80], l_urban: [0, 20], q90: [0.05, 2], q_nu: [0.5, 3], x_br: [0, 5000] };

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

/* ---------- surface model (tiled, loaded on demand) ---------- */
const MAX_TILES = 110;   // ~2.6 MB each; tiles not needed by the current task are evicted beyond this
class Model {
  constructor(meta, base) {
    Object.assign(this, { meta, base });
    this.W = meta.width; this.H = meta.height; this.x0 = meta.x0; this.y0 = meta.y0; this.res = meta.res; this.zone = meta.utm_zone;
    this.TS = meta.tile; this.GF = meta.ground_factor; this.GT = meta.tile / meta.ground_factor + 2;
    this.tiles = new Map(); this.clock = 0; this.walls = new Map(); this.streets = new Map();
    const c = [[0, this.H], [this.W, this.H], [this.W, 0], [0, 0]].map(([cc, rr]) => this.toLL(this.x0 + cc * this.res, this.y0 - rr * this.res));
    this.corners = c.map(([lon, lat]) => [lat, lon]);   // SW, SE, NE, NW
    this.bounds = [[Math.min(...c.map(p => p[1])), Math.min(...c.map(p => p[0]))], [Math.max(...c.map(p => p[1])), Math.max(...c.map(p => p[0]))]];
  }
  toXY(lon, lat) { return UTM.fwd(lon, lat, this.zone); }
  toLL(x, y) { return UTM.inv(x, y, this.zone); }
  col(x) { return (x - this.x0) / this.res; }
  row(y) { return (this.y0 - y) / this.res; }
  inside(x, y) {
    const r = Math.floor(this.row(y)), c = Math.floor(this.col(x));
    return r >= 0 && r < this.H && c >= 0 && c < this.W && (`${Math.floor(c / this.TS)}_${Math.floor(r / this.TS)}` in this.meta.tiles);
  }
  tileAt(r, c) {
    const k = `${Math.floor(c / this.TS)}_${Math.floor(r / this.TS)}`, t = this.tiles.get(k);
    if (!t) throw new Error('Map data for this area is not loaded (tile ' + k + ')');
    return t;
  }
  rc(x, y) {
    let r = Math.floor(this.row(y)), c = Math.floor(this.col(x));
    return [r < 0 ? 0 : r >= this.H ? this.H - 1 : r, c < 0 ? 0 : c >= this.W ? this.W - 1 : c];
  }
  surfRC(r, c) { return this.tileAt(r, c).S[(r % this.TS) * this.TS + (c % this.TS)] * 0.1; }
  taggedRC(r, c) { const i = (r % this.TS) * this.TS + (c % this.TS); return (this.tileAt(r, c).T[i >> 3] >> (i & 7)) & 1; }
  surf(x, y) { const [r, c] = this.rc(x, y); return this.surfRC(r, c); }
  tagged(x, y) { const [r, c] = this.rc(x, y); return this.taggedRC(r, c); }
  ground(x, y) {   // bilinear on the coarse grid; each tile carries a one-cell border
    const [r, c] = this.rc(x, y), t = this.tileAt(r, c);
    const lr = this.row(y) - (r - r % this.TS), lc = this.col(x) - (c - c % this.TS);
    let u = lc / this.GF + 0.5, v = lr / this.GF + 0.5;
    const G = t.G, n = this.GT;
    u = Math.min(Math.max(u, 0), n - 1.001); v = Math.min(Math.max(v, 0), n - 1.001);
    const i = Math.floor(u), j = Math.floor(v), fu = u - i, fv = v - j;
    const a = G[j * n + i], b = G[j * n + i + 1], cc = G[(j + 1) * n + i], d = G[(j + 1) * n + i + 1];
    return 0.1 * ((a * (1 - fu) + b * fu) * (1 - fv) + (cc * (1 - fu) + d * fu) * fv);
  }
  /* make sure every tile touching the UTM box [x0,y0,x1,y1] is loaded */
  async ensure(box, onProgress) {
    const c0 = Math.max(Math.floor(this.col(box[0]) / this.TS), 0), c1 = Math.min(Math.floor(this.col(box[2]) / this.TS), this.W / this.TS - 1);
    const r0 = Math.max(Math.floor(this.row(box[3]) / this.TS), 0), r1 = Math.min(Math.floor(this.row(box[1]) / this.TS), this.H / this.TS - 1);
    const need = [];
    for (let ty = r0; ty <= r1; ty++) for (let tx = c0; tx <= c1; tx++) { const k = `${tx}_${ty}`; if (k in this.meta.tiles) need.push(k); }
    const stamp = ++this.clock, missing = need.filter(k => !this.tiles.has(k));
    need.forEach(k => { const t = this.tiles.get(k); if (t) t.used = stamp; });
    if (missing.length) {
      const total = missing.reduce((s, k) => s + this.meta.tiles[k], 0); let got = 0;
      const queue = missing.slice();
      const worker = async () => {
        while (queue.length) {
          const k = queue.shift();
          const buf = await gunzipMaybe(await fetchBuf(this.base + `t_${k}.bin.gz`, n => { got += n; onProgress && onProgress(Math.min(got / total, 1)); }));
          this.tiles.set(k, Object.assign(decodeTile(buf, this.TS, this.GT), { used: stamp }));
        }
      };
      await Promise.all(Array.from({ length: Math.min(6, queue.length) }, worker));
    }
    if (this.tiles.size > MAX_TILES) {
      const old = [...this.tiles.entries()].filter(([, t]) => t.used !== stamp).sort((p, q) => p[1].used - q[1].used);
      for (const [k] of old) { if (this.tiles.size <= MAX_TILES) break; this.tiles.delete(k); }
    }
    return missing.length;
  }
  /* building walls (for reflections): every tile whose walls could matter for the box */
  async ensureWalls(box, onProgress) {
    if (!this.meta.walls) return 0;
    const c0 = Math.max(Math.floor(this.col(box[0]) / this.TS) - 1, 0), c1 = Math.min(Math.floor(this.col(box[2]) / this.TS) + 1, this.W / this.TS - 1);
    const r0 = Math.max(Math.floor(this.row(box[3]) / this.TS) - 1, 0), r1 = Math.min(Math.floor(this.row(box[1]) / this.TS) + 1, this.H / this.TS - 1);
    const missing = [];
    for (let ty = r0; ty <= r1; ty++) for (let tx = c0; tx <= c1; tx++) { const k = `${tx}_${ty}`; if (k in this.meta.walls && !this.walls.has(k)) missing.push(k); }
    const total = missing.reduce((s, k) => s + this.meta.walls[k], 0) || 1; let got = 0;
    await Promise.all(missing.map(async k => {
      const buf = await gunzipMaybe(await fetchBuf(this.base + `w_${k}.bin.gz`, n => { got += n; onProgress && onProgress(Math.min(got / total, 1)); }));
      const [tx, ty] = k.split('_').map(Number);
      this.walls.set(k, decodeWalls(buf, this.x0 + tx * this.TS * this.res, this.y0 - (ty + 1) * this.TS * this.res));
    }));
    return missing.length;
  }
  async ensureStreets(box, onProgress) {
    if (!this.meta.streets) return 0;
    const c0 = Math.max(Math.floor(this.col(box[0]) / this.TS) - 1, 0), c1 = Math.min(Math.floor(this.col(box[2]) / this.TS) + 1, this.W / this.TS - 1);
    const r0 = Math.max(Math.floor(this.row(box[3]) / this.TS) - 1, 0), r1 = Math.min(Math.floor(this.row(box[1]) / this.TS) + 1, this.H / this.TS - 1);
    const missing = [];
    for (let ty = r0; ty <= r1; ty++) for (let tx = c0; tx <= c1; tx++) { const k = `${tx}_${ty}`; if (k in this.meta.streets && !this.streets.has(k)) missing.push(k); }
    const total = missing.reduce((s, k) => s + this.meta.streets[k], 0) || 1; let got = 0;
    await Promise.all(missing.map(async k => {
      const buf = await gunzipMaybe(await fetchBuf(this.base + `s_${k}.bin.gz`, n => { got += n; onProgress && onProgress(Math.min(got / total, 1)); }));
      const [tx, ty] = k.split('_').map(Number);
      this.streets.set(k, decodeStreets(buf, this.x0 + tx * this.TS * this.res, this.y0 - (ty + 1) * this.TS * this.res));
    }));
    return missing.length;
  }
  wallTilesIn(box) {
    const out = [];
    for (const [k, w] of this.walls) if (!(w.box[2] < box[0] || w.box[0] > box[2] || w.box[3] < box[1] || w.box[1] > box[3])) out.push(w);
    return out;
  }
  boxLL(list, padM) {   // UTM box around lat/lon points, padded by padM metres
    const xy = list.map(p => this.toXY(+p.lon, +p.lat));
    return [Math.min(...xy.map(q => q[0])) - padM, Math.min(...xy.map(q => q[1])) - padM,
            Math.max(...xy.map(q => q[0])) + padM, Math.max(...xy.map(q => q[1])) + padM];
  }
}
/* walls: Float64Array of segments [x1, y1, x2, y2] in UTM, outside to the right (CCW rings) */
function decodeWalls(buf, ox, oy) {
  const a = new Int16Array(buf), segs = [];
  let i = 0, bx0 = Infinity, by0 = Infinity, bx1 = -Infinity, by1 = -Infinity;
  while (i < a.length) {
    const n = a[i++]; let x = a[i++], y = a[i++];
    const xs = [x], ys = [y];
    for (let k = 1; k < n; k++) { x += a[i++]; y += a[i++]; xs.push(x); ys.push(y); }
    for (let k = 0; k < n; k++) {
      const x1 = ox + xs[k] / 2, y1 = oy + ys[k] / 2, x2 = ox + xs[(k + 1) % n] / 2, y2 = oy + ys[(k + 1) % n] / 2;
      segs.push(x1, y1, x2, y2);
      bx0 = Math.min(bx0, x1); by0 = Math.min(by0, y1); bx1 = Math.max(bx1, x1); by1 = Math.max(by1, y1);
    }
  }
  return { segs: new Float64Array(segs), box: [bx0, by0, bx1, by1] };
}
/* streets: list of polylines (Float64Array [x0, y0, x1, y1, ...] in UTM); pieces meet at identical ends */
function decodeStreets(buf, ox, oy) {
  const a = new Int16Array(buf), lines = [];
  let i = 0;
  while (i < a.length) {
    const n = a[i++]; let x = a[i++], y = a[i++];
    const p = new Float64Array(2 * n); p[0] = ox + x / 2; p[1] = oy + y / 2;
    for (let k = 1; k < n; k++) { x += a[i++]; y += a[i++]; p[2 * k] = ox + x / 2; p[2 * k + 1] = oy + y / 2; }
    lines.push(p);
  }
  return lines;
}
function decodeTile(buf, TS, GT) {
  const nS = TS * TS, nG = GT * GT;
  if (buf.byteLength !== 2 * nS + 2 * nG + nS / 8) throw new Error('map tile has unexpected size');
  return { S: undelta(buf.slice(0, 2 * nS), TS, TS), G: undelta(buf.slice(2 * nS, 2 * nS + 2 * nG), GT, GT),
           T: new Uint8Array(buf.slice(2 * nS + 2 * nG)) };
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
  // split at the main edge only if it actually blocks the line of sight; an edge just below the
  // line only costs its own (small) loss - splitting there would measure other obstacles against
  // lowered sub-paths and overstate the loss of a nearly clear path
  if (best <= 0) return J(best);
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

/* ---------- single-bounce wall reflections ----------
 * Every outer building wall is a flat mirror. For a pair a-b we mirror b in each wall's line,
 * keep walls where the bounce point lies on the wall (inside its length, above ground, below
 * the roof edge read from the surface model) and both points are in front of it.
 * Path loss = free space over d1+d2 + diffraction on each leg (Deygout, as for direct links)
 * + reflection loss: Fresnel coefficient for a vertical antenna on a dielectric wall
 * (relative permittivity wall_er, 5 = concrete/brick) + a size penalty when the wall is smaller
 * than the first Fresnel zone on it + wall_extra dB for windows, balconies and roughness.
 */
function wallLoss(cosI, er) {
  const s2 = 1 - cosI * cosI, q = Math.sqrt(er - s2), g = Math.abs((cosI - q) / (cosI + q));
  return -20 * Math.log10(Math.max(g, 1e-3));
}
function reflectionPaths(M, a, b, P, keep = 5) {
  const lam = 299.792458 / P.freq, D = Math.hypot(b.x - a.x, b.y - a.y);
  const excess = Math.min(Math.max(0.3 * D, 1000), 3000);
  const semiMin = Math.sqrt(((D + excess) / 2) ** 2 - (D / 2) ** 2);
  const box = [Math.min(a.x, b.x) - semiMin, Math.min(a.y, b.y) - semiMin, Math.max(a.x, b.x) + semiMin, Math.max(a.y, b.y) + semiMin];
  const fs = d => 20 * Math.log10(Math.max(d, 1) / 1000) + 20 * Math.log10(P.freq) + 32.44;
  const cands = [];
  for (const w of M.wallTilesIn(box)) {
    const S = w.segs;
    for (let i = 0; i < S.length; i += 4) {
      const x1 = S[i], y1 = S[i + 1], ex = S[i + 2] - x1, ey = S[i + 3] - y1, L = Math.hypot(ex, ey);
      if (L < 6) continue;
      const nx = ey / L, ny = -ex / L;                       // outward normal (CCW ring)
      const sa = (a.x - x1) * nx + (a.y - y1) * ny, sb = (b.x - x1) * nx + (b.y - y1) * ny;
      if (sa < 2 || sb < 2) continue;                        // both must be in front of the wall
      const t = sa / (sa + sb);
      const rx = a.x + (b.x - 2 * sb * nx - a.x) * t, ry = a.y + (b.y - 2 * sb * ny - a.y) * t;
      const u = ((rx - x1) * ex + (ry - y1) * ey) / (L * L);
      if (u < 0.05 || u > 0.95) continue;
      const d1 = Math.hypot(rx - a.x, ry - a.y), d2 = Math.hypot(b.x - rx, b.y - ry);
      if (d1 < 15 || d2 < 15 || d1 + d2 > D + excess) continue;
      const cosI = Math.min(sa / d1, 1);
      // Fresnel zone on the wall: radius r, stretched horizontally by 1/cos(incidence)
      const r = Math.sqrt(lam * d1 * d2 / (d1 + d2));
      const wFrac = Math.min(1, L * Math.min(u, 1 - u) * 2 * cosI / r);
      const quick = fs(d1 + d2) + wallLoss(cosI, P.wall_er) - 20 * Math.log10(Math.max(wFrac, 0.05));
      cands.push({ rx, ry, nx, ny, d1, d2, cosI, r, wFrac, quick, L });
    }
  }
  cands.sort((p, q) => p.quick - q.quick);
  const found = [];
  for (const c of cands.slice(0, 250)) {
    const zR = a.z + (b.z - a.z) * c.d1 / (c.d1 + c.d2);
    const g = M.ground(c.rx, c.ry), top = M.surf(c.rx - c.nx * 2.5, c.ry - c.ny * 2.5);
    if (zR < g + 1 || zR > top - 0.5) continue;              // bounce point must be on the wall face
    const hFrac = Math.min(1, Math.min(zR - g, top - zR) * 2 / c.r);
    const R = { x: c.rx + c.nx, y: c.ry + c.ny, z: zR };
    const l1 = link(M, a, R, P.freq).diff, l2 = link(M, R, b, P.freq).diff;
    const sizePen = Math.min(20, -20 * Math.log10(Math.max(c.wFrac * hFrac, 1e-2)));
    const refl = wallLoss(c.cosI, P.wall_er) + sizePen + P.wall_extra;
    const loss = fs(c.d1 + c.d2) + l1 + l2 + refl;
    found.push({ loss, rx: c.rx, ry: c.ry, d1: c.d1, d2: c.d2, refl, diff: l1 + l2,
      angle: Math.round(Math.acos(c.cosI) * 180 / Math.PI), wall_m: Math.round(c.L), height: Math.round(zR - g) });
  }
  found.sort((p, q) => p.loss - q.loss);
  return { paths: found.slice(0, keep), candidates: cands.length, checked: Math.min(cands.length, 250) };
}

/* ---------- street-level multipath, statistical (ITU-R P.1411-12, section 4.3.1) ----------
 * Site-general model for terminals near street level, built from measurements in cities:
 *   L_NLoS(median) = 9.5 + 45 log10 f[MHz] + 40 log10(d[m]/1000) + L_urban
 *   L_urban = 0 (suburban), 6.8 (urban), 2.3 (dense urban / high-rise); location spread sigma = 7 dB.
 * ITU validity: 300-3000 MHz, d <= 3 km, antennas 1.9-3 m. Used here only if both antennas are
 * <= 10 m above ground, the path crosses built-up area (>= 5 % of samples inside buildings)
 * and d <= 5 km (3-5 km flagged as extrapolated). Returns null when it does not apply.
 */
function streetLoss(M, a, b, P) {
  if (P.freq < 300 || P.freq > 3000) return null;
  const D = Math.hypot(b.x - a.x, b.y - a.y);
  if (D < 50 || D > 5000) return null;
  const ha = a.z - M.ground(a.x, a.y), hb = b.z - M.ground(b.x, b.y);
  if (ha > 10 || hb > 10) return null;
  let inB = 0, n = 0;
  for (let d = 10; d < D - 10; d += 20) { const t = d / D, x = a.x + (b.x - a.x) * t, y = a.y + (b.y - a.y) * t; n++; if (M.surf(x, y) - M.ground(x, y) > 2) inB++; }
  if (!n || inB / n < 0.05) return null;
  const fspl = 20 * Math.log10(D / 1000) + 20 * Math.log10(P.freq) + 32.44;
  const L = 9.5 + 45 * Math.log10(P.freq) + 40 * Math.log10(D / 1000) + P.l_urban;
  return { loss: Math.max(L, fspl), built_share: Math.round(inB / n * 100), extrapolated: D > 3000 };
}

/* ---------- street corridors: Berg's recursive street model ----------
 * The signal runs along straight street segments; each corner adds loss that grows with the
 * distances on both sides of it ("illusory distance"), as in the micro-cell model of
 * ETSI TR 125 942 (3GPP TR 25.942) after J.-E. Berg, PIMRC 1995:
 *   k_j = k_(j-1) + d_(j-1) * q(theta_j),  d_j = k_j * s_j + d_(j-1),  k_0 = 1, d_0 = 0
 *   q(theta) = q90 * (theta / 90)^nu  (q = q90 = 0.5 for a 90-degree corner, as in ETSI)
 *   L = 20 log10(4 pi d_n / lambda) + 20 log10(prod D(s_j)),  D(s) = s / x_br beyond x_br (300 m)
 *       (slope 2 before the 300 m breakpoint and 4 after it, applied per straight segment)
 * A straight segment is a line-of-sight run along the street: the route stays in it while it keeps
 * within +-10 m (about half a street) of the line from the segment's start; jogs at junctions,
 * medians and kinks don't count as corners. A corner's angle is measured between the segment's
 * direction and the new direction; corners under 15 degrees cost nothing. The route is searched only in a
 * corridor around the straight line between the points. Each point joins the street network at
 * a street spot within 200 m that it can see (or, failing that, with the diffraction loss of
 * that short leg added); the leg itself is the first / last segment.
 */
const STRAIGHT_DEG = 15, MAX_TURNS = 6, HALF_STREET = 10;
function bergRoute(M, a, b, P) {
  const lam = 299.792458 / P.freq, D = Math.hypot(b.x - a.x, b.y - a.y);
  if (D < 30) return { none: true, reason: 'points closer than 30 m' };
  const hw = Math.min(Math.max(0.25 * D, 300), 1500), ATT = 200;
  const ux = (b.x - a.x) / D, uy = (b.y - a.y) / D;
  const inCorr = (x, y) => { const t = (x - a.x) * ux + (y - a.y) * uy, n = Math.abs(-(x - a.x) * uy + (y - a.y) * ux);
    return t > -ATT && t < D + ATT && n < hw; };
  // graph of densified street vertices inside the corridor
  const id = new Map(), X = [], Y = [], adj = [];
  const node = (x, y) => { const k = Math.round(x * 2) * 4e6 + Math.round(y * 2); let i = id.get(k);
    if (i === undefined) { i = X.length; id.set(k, i); X.push(x); Y.push(y); adj.push([]); } return i; };
  const box = [Math.min(a.x, b.x) - hw - ATT, Math.min(a.y, b.y) - hw - ATT, Math.max(a.x, b.x) + hw + ATT, Math.max(a.y, b.y) + hw + ATT];
  for (const lines of M.streets.values()) for (const p of lines) {
    let anyIn = false;
    for (let k = 0; k < p.length; k += 2) if (inCorr(p[k], p[k + 1])) { anyIn = true; break; }
    if (!anyIn) continue;
    let prev = node(p[0], p[1]);
    for (let k = 2; k < p.length; k += 2) {
      const x0 = p[k - 2], y0 = p[k - 1], x1 = p[k], y1 = p[k + 1], len = Math.hypot(x1 - x0, y1 - y0);
      const m = Math.max(1, Math.ceil(len / 25));          // vertices every <= 25 m, for joining points
      for (let s = 1; s <= m; s++) {
        const cur = s === m ? node(x1, y1) : node(x0 + (x1 - x0) * s / m, y0 + (y1 - y0) * s / m);
        if (cur !== prev) { adj[prev].push(cur); adj[cur].push(prev); }
        prev = cur;
      }
    }
  }
  if (!X.length) return { none: true, reason: 'no streets in the corridor' };
  const q = th => th < STRAIGHT_DEG ? 0 : P.q90 * Math.pow(th / 90, P.q_nu);
  const Dbp = s => P.x_br > 0 && s > P.x_br ? s / P.x_br : 1;   // x_br = 0 switches the breakpoint off
  const turn = (dx1, dy1, dx2, dy2) => { const c = (dx1 * dx2 + dy1 * dy2) / (Math.hypot(dx1, dy1) * Math.hypot(dx2, dy2) || 1);
    return Math.acos(Math.max(-1, Math.min(1, c))) * 180 / Math.PI; };
  // joining legs: street spots within ATT m that the point can see (or the least obstructed ones)
  const legs = (pt) => {
    const out = [];
    for (let i = 0; i < X.length; i++) {
      const d = Math.hypot(X[i] - pt.x, Y[i] - pt.y);
      if (d > ATT) continue;
      out.push([i, d]);
    }
    out.sort((p, r) => p[1] - r[1]);
    const res = [];
    for (const [i, d] of out.slice(0, 60)) {
      const S = { x: X[i], y: Y[i], z: M.ground(X[i], Y[i]) + 1.5 };
      const pen = d < 12 ? 0 : link(M, pt, S, P.freq).diff;
      if (pen <= 25) res.push({ i, d, pen });
    }
    return res;
  };
  // stays straight if w is within HALF_STREET of the line from the segment start through v (or the segment is still short)
  const offLine = (sx, sy, vx, vy, wx, wy) => {
    const lx = vx - sx, ly = vy - sy, l = Math.hypot(lx, ly);
    if (l < HALF_STREET) return 0;
    const t = ((wx - sx) * lx + (wy - sy) * ly) / l;
    if (t < l - 1) return Infinity;                       // going backwards = a corner
    return Math.abs(-(wx - sx) * ly + (wy - sy) * lx) / l;
  };
  const startLegs = legs(a), endLegs = legs(b);
  if (!startLegs.length || !endLegs.length)
    return { none: true, reason: `${!startLegs.length ? 'first' : 'second'} point is not in view of a street within ${ATT} m` };
  const exitAt = new Map(endLegs.map(e => [e.i, e]));
  // Dijkstra over (node, came-from node); label = illusory-distance state
  const heap = [], best = new Map(), from = new Map();
  const push = (st) => { heap.push(st); let i = heap.length - 1;
    while (i > 0) { const p = (i - 1) >> 1; if (heap[p].pr <= heap[i].pr) break; [heap[p], heap[i]] = [heap[i], heap[p]]; i = p; } };
  const pop = () => { const top = heap[0], last = heap.pop();
    if (heap.length) { heap[0] = last; let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i;
      if (l < heap.length && heap[l].pr < heap[m].pr) m = l; if (r < heap.length && heap[r].pr < heap[m].pr) m = r;
      if (m === i) break; [heap[m], heap[i]] = [heap[i], heap[m]]; i = m; } } return top; };
  const prio = st => (st.k * st.s + st.d) * st.Dp * Dbp(st.s) * Math.pow(10, st.pen / 20);
  for (const L0 of startLegs) {
    const st = { node: L0.i, prevX: a.x, prevY: a.y, sx: a.x, sy: a.y, k: 1, d: 0, s: L0.d, Dp: 1, pen: L0.pen, turns: 0, key: L0.i + ':a' };
    st.pr = prio(st); push(st); from.set(st.key, null);
  }
  let bestEnd = null, pops = 0;
  while (heap.length && pops < 400000) {
    const st = pop(); pops++;
    if (bestEnd && st.pr >= bestEnd.pr) break;
    if (best.has(st.key) && best.get(st.key) <= st.pr) continue;
    best.set(st.key, st.pr);
    const e = exitAt.get(st.node);
    if (e) {   // try finishing here: last leg to b
      const v0 = st.node, straight = e.d < 1 || offLine(st.sx, st.sy, X[v0], Y[v0], b.x, b.y) <= HALF_STREET;
      const th = straight ? 0 : turn(X[v0] - st.sx, Y[v0] - st.sy, b.x - X[v0], b.y - Y[v0]);
      let k = st.k, d = st.d, s = st.s, Dp = st.Dp, turns = st.turns;
      if (q(th) > 0) { if (turns >= (P.max_turns || MAX_TURNS)) continue; d = k * s + d; Dp *= Dbp(s); k = k + d * q(th); s = e.d; turns++; } else s += e.d;
      const fin = { d: k * s + d, Dp: Dp * Dbp(s), pen: st.pen + e.pen, turns, key: st.key, node: st.node };
      fin.pr = fin.d * fin.Dp * Math.pow(10, fin.pen / 20);
      if (!bestEnd || fin.pr < bestEnd.pr) bestEnd = fin;
    }
    const v = st.node;
    for (const w of adj[v]) {
      if (!inCorr(X[w], Y[w]) || w === st.prevNode) continue;
      const len = Math.hypot(X[w] - X[v], Y[w] - Y[v]);
      let k = st.k, d = st.d, s = st.s, Dp = st.Dp, turns = st.turns, sx = st.sx, sy = st.sy;
      if (offLine(sx, sy, X[v], Y[v], X[w], Y[w]) > HALF_STREET) {
        const th = turn(X[v] - sx, Y[v] - sy, X[w] - X[v], Y[w] - Y[v]);
        if (q(th) > 0) { if (turns >= (P.max_turns || MAX_TURNS)) continue; d = k * s + d; Dp *= Dbp(s); k = k + d * q(th); s = 0; turns++; }
        sx = X[v]; sy = Y[v];
      }
      s += len;
      const nst = { node: w, prevNode: v, prevX: X[v], prevY: Y[v], sx, sy, k, d, s, Dp, pen: st.pen, turns, key: w + ':' + Math.round(sx) + ',' + Math.round(sy) };
      nst.pr = prio(nst);
      if (best.has(nst.key) && best.get(nst.key) <= nst.pr) continue;
      from.set(nst.key, st.key); push(nst);
    }
  }
  if (!bestEnd) return { none: true, reason: `no street route with at most ${P.max_turns || MAX_TURNS} corners inside the corridor` };
  // rebuild the route
  const route = [[b.x, b.y]];
  for (let k = bestEnd.key; k; k = from.get(k)) { const n = +k.split(':')[0]; route.push([X[n], Y[n]]); }
  route.push([a.x, a.y]); route.reverse();
  let len = 0; for (let i = 1; i < route.length; i++) len += Math.hypot(route[i][0] - route[i - 1][0], route[i][1] - route[i - 1][1]);
  const loss = 20 * Math.log10(4 * Math.PI * bestEnd.d * bestEnd.Dp / lam) + bestEnd.pen;
  // how much of the route runs between building walls (both sides within 30 m)
  let walled = 0, nS = 0;
  for (let i = 1; i < route.length; i++) {
    const [x0, y0] = route[i - 1], [x1, y1] = route[i], l = Math.hypot(x1 - x0, y1 - y0); if (l < 1) continue;
    const nx = -(y1 - y0) / l, ny = (x1 - x0) / l;
    for (let t = 5; t < l; t += 20) {
      const x = x0 + (x1 - x0) * t / l, y = y0 + (y1 - y0) * t / l; nS++;
      const side = sg => { for (let o = 6; o <= 30; o += 4) { const px = x + sg * nx * o, py = y + sg * ny * o; if (M.inside(px, py) && M.surf(px, py) - M.ground(px, py) > 2) return true; } return false; };
      if (side(1) && side(-1)) walled++;
    }
  }
  return { loss, illusory_m: bestEnd.d * bestEnd.Dp, turns: bestEnd.turns, route_m: len, join_db: bestEnd.pen,
    walled: nS ? Math.round(walled / nS * 100) : 0, route: route.map(([x, y]) => { const [lon, lat] = M.toLL(x, y); return [Math.round(lat * 1e6) / 1e6, Math.round(lon * 1e6) / 1e6]; }) };
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

function computeLinks(M, pts, P, withRefl, onProgress, withCorr) {
  const budget = P.tx_dbm + P.gain_tx + P.gain_rx - P.losses, rows = [], total = pts.length * (pts.length - 1) / 2;
  for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++) {
    const L = link(M, pts[i], pts[j], P.freq);
    let loss = L.loss, refl = null;
    if (withRefl) {
      const R = reflectionPaths(M, pts[i], pts[j], P);
      // incoherent power sum of the direct path and the best reflected paths
      const pw = [L.loss, ...R.paths.map(p => p.loss)].reduce((s, l) => s + Math.pow(10, -l / 10), 0);
      const det = -10 * Math.log10(pw), st = streetLoss(M, pts[i], pts[j], P);
      const br = withCorr ? bergRoute(M, pts[i], pts[j], P) : null;
      loss = det;
      if (st && st.loss < loss) loss = st.loss;
      // a route without corners is plain line of sight along one street - the direct estimate covers it
      if (br && !br.none && br.turns >= 1 && br.loss < loss) loss = br.loss;
      refl = { found: R.paths.length, candidates: R.candidates, direct_margin: r1(budget - L.loss - P.sens_dbm),
        gain_db: r1(L.loss - loss), walls_margin: r1(budget - det - P.sens_dbm), walls_gain: r1(L.loss - det),
        street: st ? { margin: r1(budget - st.loss - P.sens_dbm), built_share: st.built_share, extrapolated: st.extrapolated } : null,
        corridor: !withCorr ? null : br.none ? { none: true, reason: br.reason }
          : { margin: r1(budget - br.loss - P.sens_dbm), turns: br.turns, route_m: Math.round(br.route_m), walled: br.walled, join_db: r1(br.join_db), route: br.route },
        used: br && !br.none && br.turns >= 1 && br.loss === loss ? 'corridor' : st && st.loss === loss ? 'street' : (L.loss - det >= 0.5 ? 'walls' : 'direct'),
        best: R.paths[0] ? (() => { const p = R.paths[0], [lon, lat] = M.toLL(p.rx, p.ry);
          return { lat: Math.round(lat * 1e6) / 1e6, lon: Math.round(lon * 1e6) / 1e6, margin: r1(budget - p.loss - P.sens_dbm),
            extra_m: Math.round(p.d1 + p.d2 - L.dist_km * 1000), refl_loss: r1(p.refl), diff: r1(p.diff), angle: p.angle, wall_m: p.wall_m, height: p.height }; })() : null,
        others: R.paths.slice(1).map(p => { const [lon, lat] = M.toLL(p.rx, p.ry); return { lat, lon, margin: r1(budget - p.loss - P.sens_dbm) }; }) };
      if (onProgress) onProgress(rows.length / total);
    }
    const rx = budget - loss, margin = rx - P.sens_dbm;
    rows.push({ a: pts[i].name, b: pts[j].name, dist_km: Math.round(L.dist_km * 100) / 100, los: L.los,
      fresnel: Math.round(L.fresnel * 100) / 100, fspl: r1(L.fspl), diff: r1(L.diff),
      model: !refl || refl.used === 'direct' ? L.model : refl.used === 'walls' ? L.model + ' + wall reflection'
        : refl.used === 'corridor' ? "street corridor (Berg's recursive model)" : 'street multipath (ITU-R P.1411)' + (refl.street.extrapolated ? ', extrapolated beyond 3 km' : ''),
      rx_dbm: r1(rx), margin: r1(margin), status: classify(margin), refl,
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
    const x = M.x0 + (c0 + c + 0.5) * M.res, y = M.y0 - (r0 + r + 0.5) * M.res, s = M.surfRC(r0 + r, c0 + c);
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
    const tg = M.taggedRC(r0 + r, c0 + c);
    if ((hb[i] >= P.min_roof && tg) || hb[i] < 2) cands.push([S[i], r, c]);
  }
  cands.sort((p, q) => q[0] - p[0]);
  const up = P.tx_dbm + P.gain_tx + P.rep_gain - P.losses, dn = P.rep_tx_dbm + P.rep_gain + P.gain_rx - P.losses;
  const budget = Math.min(up, dn), res = [], list = cands.slice(0, P.max_candidates);
  list.forEach(([s, r, c], k) => {
    const x = M.x0 + (c0 + c + 0.5) * M.res, y = M.y0 - (r0 + r + 0.5) * M.res;
    const surfC = M.surfRC(r0 + r, c0 + c), bh = hb[r * w + c];
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
  const links = pts.length >= 2 ? computeLinks(M, pts, P, !!body.reflections, p => onProgress && onProgress(p, 'reflections'), !!(body.reflections && body.corridors)) : [];
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

/* ---------- viewshed ("visibility cloak"), one or many observers ----------
 * All observers share one lat/lon grid (cell ~cellM metres). From each observer, rays go to every
 * edge cell of its own square; along a ray the terrain+buildings surface is sampled every 5 m.
 * A cell counts as visible when a target standing `targetH` metres above the surface there
 * (ground or roof) is above the highest elevation angle met so far. 4/3 earth curvature.
 * Rays are straight in lat/lon; their UTM positions come from the exact projection every 200 m
 * (linear in between), so cell position and sampled terrain agree to well under a metre.
 * Result: bits[i] has bit k set if observer k sees ALL of cell i (every sample in it); cov[i] bit k if cell i is inside
 * observer k's radius and inside the model.
 */
const mLatAt = lat => 111132.92 - 559.82 * Math.cos(2 * lat * Math.PI / 180);
const mLonAt = lat => 111412.84 * Math.cos(lat * Math.PI / 180);

function viewshed(M, body, onProgress) {
  const raw = body.points || (body.point ? [body.point] : []);
  if (!raw.length) throw new Error('No point selected');
  if (raw.length > 24) throw new Error('At most 24 points at once');
  const radius = Math.min(Math.max(+body.radius_m || 10000, 500), 30000);
  const targetH = Math.min(Math.max(+body.target_h || 0, 0), 300);
  const obs = placePoints(M, raw);
  // common grid
  let south = Infinity, north = -Infinity, west = Infinity, east = -Infinity;
  for (const o of obs) {
    const dla = radius / mLatAt(o.lat), dlo = radius / mLonAt(o.lat + Math.sign(o.lat) * dla);
    south = Math.min(south, o.lat - dla); north = Math.max(north, o.lat + dla);
    west = Math.min(west, o.lon - dlo); east = Math.max(east, o.lon + dlo);
  }
  const latc = (south + north) / 2;
  const extent = Math.max((north - south) * mLatAt(latc), (east - west) * mLonAt(latc));
  const cellM = Math.max(+body.cell_m || 0, 10, Math.ceil(extent / 2200));
  const dLat = cellM / mLatAt(latc), dLon = cellM / mLonAt(latc);
  const H = Math.ceil((north - south) / dLat), W = Math.ceil((east - west) / dLon);
  north = south + H * dLat; east = west + W * dLon;
  const bits = new Uint32Array(W * H), cov = new Uint32Array(W * H), blk = new Uint32Array(W * H);
  const STEP_V = cellM <= 20 ? 5 : 10, KNOT = 200;   // 10 m = model resolution; finer only for small areas
  const stats = obs.map(() => ({ vis: 0, tot: 0 }));

  obs.forEach((o, k) => {
    const bit = 1 << k;
    const r0 = Math.floor((north - o.lat) / dLat), c0 = Math.floor((o.lon - west) / dLon);
    const hr = Math.ceil(radius / cellM) + 1;
    const edge = [];
    for (let t = -hr; t <= hr; t++) edge.push([r0 - hr, c0 + t], [r0 + hr, c0 + t], [r0 + t, c0 - hr], [r0 + t, c0 + hr]);
    for (const [er, ec] of edge) {
      const tlat = north - (er + 0.5) * dLat, tlon = west + (ec + 0.5) * dLon;
      const [tx, ty] = M.toXY(tlon, tlat);
      const len = Math.hypot(tx - o.x, ty - o.y);
      const nK = Math.max(Math.ceil(len / KNOT), 1);
      // exact UTM at knots along the lat/lon-straight ray
      const kx = new Float64Array(nK + 1), ky = new Float64Array(nK + 1);
      for (let j = 0; j <= nK; j++) { const f = j / nK; const [x, y] = M.toXY(o.lon + (tlon - o.lon) * f, o.lat + (tlat - o.lat) * f); kx[j] = x; ky[j] = y; }
      let maxSlope = -Infinity;
      for (let d = STEP_V; d <= Math.min(len, radius); d += STEP_V) {
        const f = d / len, u = f * nK, j = Math.min(Math.floor(u), nK - 1), g = u - j;
        const x = kx[j] + (kx[j + 1] - kx[j]) * g, y = ky[j] + (ky[j + 1] - ky[j]) * g;
        if (!M.inside(x, y)) break;
        const dist = Math.hypot(x - o.x, y - o.y);
        if (dist > radius) break;
        const r = Math.floor((north - (o.lat + (tlat - o.lat) * f)) / dLat), c = Math.floor((o.lon + (tlon - o.lon) * f - west) / dLon);
        if (r < 0 || r >= H || c < 0 || c >= W) break;
        const s = M.surf(x, y) - dist * dist / (2 * R_EFF);
        const i = r * W + c;
        if (!(cov[i] & bit)) { cov[i] |= bit; stats[k].tot++; }
        if ((s + targetH - o.z) / dist >= maxSlope) bits[i] |= bit; else blk[i] |= bit;
        if (dist > 10) { const sl = (s - o.z) / dist; if (sl > maxSlope) maxSlope = sl; }
      }
    }
    if (onProgress) onProgress((k + 1) / obs.length);
  });
  // strict rule: a square counts as visible only if every sample in it was visible
  for (let i = 0; i < bits.length; i++) if (bits[i]) {
    bits[i] &= ~blk[i];
    for (let b = bits[i], k = 0; b; b >>>= 1, k++) if (b & 1) stats[k].vis++;
  }
  const a = cellM * cellM / 1e6;
  return { bits, cov, W, H, south, north, west, east, cell_m: cellM, radius_m: radius, target_h: targetH,
    observers: obs.map((o, k) => ({ name: o.name, lat: o.lat, lon: o.lon, antenna_asl: r1(o.z),
      visible_km2: r1(stats[k].vis * a), share: stats[k].tot ? Math.round(stats[k].vis / stats[k].tot * 1000) / 10 : 0 })) };
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
async function loadModel(base) {
  const res = await fetch(base + 'meta.json');
  if (!res.ok) throw new Error('Cannot load meta.json (' + res.status + ')');
  return new Model(await res.json(), base);
}

/* tiles each request needs (UTM box) */
function needBox(M, type, body) {
  const P = (() => { try { return params(body.params); } catch (e) { return DEFAULTS; } })();
  if (type === 'elevation') return M.boxLL([body], 30);
  if (type === 'profile') return M.boxLL([body.a, body.b], 50);
  if (type === 'analyze') {
    let pad = 50 + (P.repeaters ? P.search_pad : 0);
    if (body.reflections) {   // ellipses around every link (see reflectionPaths)
      const xy = (body.points || []).map(p => M.toXY(+p.lon, +p.lat));
      for (let i = 0; i < xy.length; i++) for (let j = i + 1; j < xy.length; j++) {
        const D = Math.hypot(xy[i][0] - xy[j][0], xy[i][1] - xy[j][1]), e = Math.min(Math.max(0.3 * D, 1000), 3000);
        pad = Math.max(pad, Math.sqrt(((D + e) / 2) ** 2 - (D / 2) ** 2) + 50);
        if (body.corridors) pad = Math.max(pad, Math.min(Math.max(0.25 * D, 300), 1500) + 250);
      }
    }
    return M.boxLL(body.points || [], pad);
  }
  if (type === 'viewshed') return M.boxLL(body.points || [body.point], 50 + Math.min(Math.max(+body.radius_m || 10000, 500), 30000));
  return null;
}

/* ---------- worker glue ---------- */
if (typeof self !== 'undefined' && typeof importScripts === 'function') {
  let M = null, chain = Promise.resolve();
  const handle = async ({ id, type, body }) => {
    const reply = (ok, data, transfer) => self.postMessage({ id, ok, data }, transfer || []);
    try {
      if (type === 'init') {
        M = await loadModel(body.base);
        return reply(true, { bounds: M.bounds, corners: M.corners, defaults: DEFAULTS, tiles: Object.keys(M.meta.tiles).length, total_mb: Math.round(M.meta.total_bytes / 1e5) / 10 });
      }
      if (!M) throw new Error('Model not loaded yet');
      const box = needBox(M, type, body || {});
      if (box) {
        const pts = type === 'elevation' ? [body] : type === 'profile' ? [body.a, body.b] : (body.points || [body.point]);
        for (const p of pts) { const [x, y] = M.toXY(+p.lon, +p.lat); if (!M.inside(x, y)) throw new Error(`${p.name || 'Point'} is outside the map area of the surface model`); }
        await M.ensure(box, p => self.postMessage({ type: 'progress', what: 'tiles', p }));
        if (type === 'analyze' && body.reflections) await M.ensureWalls(box, p => self.postMessage({ type: 'progress', what: 'tiles', p }));
        if (type === 'analyze' && body.reflections && body.corridors) await M.ensureStreets(box, p => self.postMessage({ type: 'progress', what: 'tiles', p }));
      }
      if (type === 'analyze') reply(true, analyze(M, body, (p, what) => self.postMessage({ type: 'progress', what: what || 'repeaters', p })));
      else if (type === 'profile') reply(true, profile(M, body));
      else if (type === 'elevation') reply(true, elevation(M, +body.lat, +body.lon));
      else if (type === 'viewshed') {
        const v = viewshed(M, body, p => self.postMessage({ type: 'progress', what: 'viewshed', p }));
        reply(true, v, [v.bits.buffer, v.cov.buffer]);
      } else throw new Error('unknown request ' + type);
    } catch (e) { reply(false, e.message || String(e)); }
  };
  // one request at a time, so tile loading and eviction never interleave
  self.onmessage = ev => { chain = chain.then(() => handle(ev.data)); };
}
if (typeof module !== 'undefined') module.exports = { UTM, Model, decodeWalls, decodeStreets, bergRoute, reflectionPaths, wallLoss, streetLoss, analyze, profile, elevation, viewshed, placePoints, params, link, loadModel, gunzipMaybe, undelta, decodeTile, needBox, DEFAULTS };
