"""Point-to-point radio visibility + repeater search over a terrain+buildings surface.

Usage:
  python3 radio_los.py points.csv --surface out/kyiv_surface.tif --ground out/kyiv_ground.tif \
      --tagged out/kyiv_height_is_tagged.tif --freq 446 --out results/

points.csv columns: name, lat, lon, height_m [, on_roof]
  height_m : antenna height in metres, above the ground (on_roof = no, default)
             or above the roof/surface at that spot (on_roof = yes)

For every pair it reports:
  - line of sight (no obstacle cuts the straight line)
  - worst Fresnel clearance (fraction of the 1st Fresnel zone left clear; >= 0.6 is "clean")
  - estimated path loss = free space + knife-edge diffraction (ITU-R P.526, Deygout, 3 edges)
  - link margin for the radio parameters below
Earth curvature uses the standard 4/3 effective radius.
"""
import argparse, csv, json, math, os
import numpy as np
import rasterio
from pyproj import Transformer

R_EFF = 6371000.0 * 4.0 / 3.0
STEP = 5.0  # metres between profile samples


def load_points(path):
    pts = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            r = {k.strip().lower(): (v or "").strip() for k, v in r.items()}
            pts.append(dict(name=r["name"], lat=float(r["lat"]), lon=float(r["lon"]),
                            h=float(r.get("height_m") or 1.5),
                            roof=r.get("on_roof", "no").lower() in ("yes", "y", "1", "true", "так")))
    return pts


def read_m(ds):
    """Read band 1 in metres (the shipped model stores int16 decimetres with scale 0.1)."""
    a = ds.read(1)
    sc, off = ds.scales[0] or 1.0, ds.offsets[0] or 0.0
    if a.dtype.kind in "iu" or sc != 1.0 or off != 0.0:
        a = a.astype(np.float32) * np.float32(sc) + np.float32(off)
    return a


class Surface:
    def __init__(self, surface, ground, tagged=None):
        self.ds = rasterio.open(surface)
        self.S = read_m(self.ds)
        self.G = read_m(rasterio.open(ground))
        self.Tg = rasterio.open(tagged).read(1) if tagged and os.path.exists(tagged) else None
        self.T = self.ds.transform
        self.to_xy = Transformer.from_crs("EPSG:4326", self.ds.crs, always_xy=True)
        self.to_ll = Transformer.from_crs(self.ds.crs, "EPSG:4326", always_xy=True)

    def rc(self, x, y):
        c = (np.asarray(x) - self.T.c) / self.T.a
        r = (np.asarray(y) - self.T.f) / self.T.e
        return r, c

    def sample(self, arr, x, y):
        r, c = self.rc(x, y)
        r = np.clip(np.floor(r).astype(int), 0, arr.shape[0] - 1)   # pixel containing the point
        c = np.clip(np.floor(c).astype(int), 0, arr.shape[1] - 1)
        return arr[r, c]

    def inside(self, x, y):
        r, c = self.rc(x, y)
        return 0 <= r < self.S.shape[0] and 0 <= c < self.S.shape[1]


def J(v):
    """ITU-R P.526 knife-edge diffraction loss (dB)."""
    v = np.asarray(v, float)
    out = np.where(v > -0.78, 6.9 + 20 * np.log10(np.sqrt((v - 0.1) ** 2 + 1) + v - 0.1), 0.0)
    return out


def deygout(d, z, h_a, h_b, lam, depth=0):
    """d: distances from start (m), z: obstacle heights incl. curvature; h_a/h_b endpoint antenna heights."""
    if len(d) < 3 or depth > 2:
        return 0.0
    D = d[-1] - d[0]
    d1 = d[1:-1] - d[0]
    d2 = D - d1
    line = h_a + (h_b - h_a) * d1 / D
    h = z[1:-1] - line
    v = h * np.sqrt(2 * D / (lam * d1 * d2))
    i = int(np.argmax(v))
    if v[i] <= -0.78:
        return 0.0
    k = i + 1
    loss = float(J(v[i]))
    loss += deygout(d[:k + 1], z[:k + 1], h_a, z[k], lam, depth + 1)
    loss += deygout(d[k:], z[k:], z[k], h_b, lam, depth + 1)
    return loss


def hata_urban(freq_mhz, hb, hm, d_km):
    """Okumura-Hata, urban, medium-size city. Valid 150-1500 MHz, hb 30-200 m, hm 1-10 m, 1-20 km."""
    hb = min(max(hb, 30.0), 200.0); hm = min(max(hm, 1.0), 10.0); d_km = max(d_km, 1.0)
    lf = math.log10(freq_mhz)
    a_hm = (1.1 * lf - 0.7) * hm - (1.56 * lf - 0.8)
    return 69.55 + 26.16 * lf - 13.82 * math.log10(hb) - a_hm + (44.9 - 6.55 * math.log10(hb)) * math.log10(d_km)


def link(srf, a, b, freq_mhz, ex_radius=10.0):
    lam = 299.792458 / freq_mhz
    D = math.hypot(b["x"] - a["x"], b["y"] - a["y"])
    n = max(int(D / STEP), 2)
    t = np.linspace(0, 1, n + 1)
    xs = a["x"] + (b["x"] - a["x"]) * t
    ys = a["y"] + (b["y"] - a["y"]) * t
    d = t * D
    z = srf.sample(srf.S, xs, ys).astype(float)
    # ignore the antenna's own building/ground spot
    z[(d < ex_radius) | (d > D - ex_radius)] = -1e4
    z += d * (D - d) / (2 * R_EFF)  # earth bulge
    ha, hb = a["z"], b["z"]
    line = ha + (hb - ha) * t
    with np.errstate(divide="ignore", invalid="ignore"):
        f1 = np.sqrt(lam * d * (D - d) / np.maximum(D, 1))
        clr = (line - z) / np.where(f1 > 0, f1, np.inf)
    mid = (d > ex_radius) & (d < D - ex_radius)
    if mid.any():
        i = int(np.argmin(np.where(mid, clr, np.inf)))
        worst = float(clr[i])
        worst_at = (xs[i], ys[i], float(line[i] - z[i]))
        tagged = bool(srf.Tg is not None and srf.sample(srf.Tg, xs[i], ys[i]) > 0)
        blocked_by_building = float(srf.sample(srf.S, xs[i], ys[i]) - srf.sample(srf.G, xs[i], ys[i])) > 2
    else:
        worst, worst_at, tagged, blocked_by_building = 99.0, None, True, False
    fspl = 20 * math.log10(max(D, 1) / 1000) + 20 * math.log10(freq_mhz) + 32.44
    zz = z.copy(); zz[~mid] = -1e4
    diff = deygout(d, zz, ha, hb, lam)
    loss, model = fspl + diff, "terrain+buildings"
    # High site <-> street-level radio: urban statistical model also applies
    # (accounts for street reflections that pure diffraction ignores).
    ga = a["z"] - float(srf.sample(srf.G, a["x"], a["y"]))
    gb = b["z"] - float(srf.sample(srf.G, b["x"], b["y"]))
    hi, lo = max(ga, gb), min(ga, gb)
    if hi >= 30 and lo <= 10 and diff > 0:
        hl = hata_urban(freq_mhz, hi, lo, D / 1000)
        if hl < loss:
            loss, model = max(hl, fspl), "urban (Hata)"
    return dict(dist_km=D / 1000, los=worst > 0, fresnel=worst, fspl=fspl, diff=diff,
                loss=loss, model=model, worst_at=worst_at, obstacle_height_known=tagged,
                obstacle_is_building=blocked_by_building)


def place(srf, p):
    x, y = srf.to_xy.transform(p["lon"], p["lat"])
    if not srf.inside(x, y):
        raise SystemExit(f"Point {p['name']} is outside the surface model area")
    base = srf.sample(srf.S if p["roof"] else srf.G, x, y)
    return dict(p, x=x, y=y, z=float(base) + p["h"], base=float(base))


def classify(margin):
    return "good" if margin >= 20 else "marginal" if margin >= 10 else "weak" if margin >= 0 else "no link"


def compute_links(srf, pts, args):
    budget = args.tx_dbm + args.gain_tx + args.gain_rx - args.losses
    rows = []
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            L = link(srf, pts[i], pts[j], args.freq)
            rx = budget - L["loss"]
            margin = rx - args.sens_dbm
            rows.append(dict(a=pts[i]["name"], b=pts[j]["name"], **L, rx_dbm=rx, margin=margin,
                             status=classify(margin)))
    return rows, budget


def run(args):
    srf = Surface(args.surface, args.ground, args.tagged)
    pts = [place(srf, p) for p in load_points(args.points)]
    os.makedirs(args.out, exist_ok=True)
    rows, budget = compute_links(srf, pts, args)
    with open(os.path.join(args.out, "links.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["from", "to", "distance_km", "line_of_sight", "fresnel_clearance", "free_space_dB",
                    "diffraction_dB", "model", "rx_level_dBm", "margin_dB", "status", "worst_obstacle",
                    "obstacle_height_from_osm_tag"])
        for r in rows:
            w.writerow([r["a"], r["b"], f"{r['dist_km']:.2f}", "yes" if r["los"] else "no",
                        f"{r['fresnel']:.2f}", f"{r['fspl']:.1f}", f"{r['diff']:.1f}", r["model"], f"{r['rx_dbm']:.1f}",
                        f"{r['margin']:.1f}", r["status"],
                        "building" if r["obstacle_is_building"] else "terrain",
                        "yes" if r["obstacle_height_known"] else "estimated"])

    reps = []
    if args.repeaters:
        reps = repeater_search(srf, pts, args, budget)
        write_repeaters_csv(reps, pts, os.path.join(args.out, "repeater_candidates.csv"))
    write_geo(srf, pts, rows, reps, args)
    return pts, rows, reps


def profile(srf, a, b, freq_mhz, n=400):
    """Height profile for plotting: distances, ground, surface (with earth bulge), line of sight, Fresnel."""
    lam = 299.792458 / freq_mhz
    D = math.hypot(b["x"] - a["x"], b["y"] - a["y"])
    t = np.linspace(0, 1, n + 1)
    xs = a["x"] + (b["x"] - a["x"]) * t
    ys = a["y"] + (b["y"] - a["y"]) * t
    d = t * D
    bulge = d * (D - d) / (2 * R_EFF)
    # max of surface over each segment so thin buildings are not missed
    seg = max(int(D / n / STEP), 1)
    tt = np.linspace(0, 1, n * seg + 1)
    Sx = srf.sample(srf.S, a["x"] + (b["x"] - a["x"]) * tt, a["y"] + (b["y"] - a["y"]) * tt)
    surf = np.array([Sx[max(k * seg - seg // 2, 0): k * seg + seg // 2 + 1].max() for k in range(n + 1)])
    ground = srf.sample(srf.G, xs, ys)
    line = a["z"] + (b["z"] - a["z"]) * t
    f1 = np.sqrt(lam * d * (D - d) / max(D, 1))
    r2 = lambda v: [round(float(x), 1) for x in v]
    return dict(d=r2(d), ground=r2(ground + bulge), surface=r2(surf + bulge), line=r2(line),
                f_lo=r2(line - f1), f_hi=r2(line + f1), za=round(a["z"], 1), zb=round(b["z"], 1))


def repeater_search(srf, pts, args, budget, top=15):
    """Candidates: rooftops of buildings with known (tagged) height >= min_roof, plus terrain highs."""
    from scipy import ndimage
    xs = [p["x"] for p in pts]; ys = [p["y"] for p in pts]
    pad = args.search_pad
    r0, c0 = srf.rc(min(xs) - pad, max(ys) + pad)
    r1, c1 = srf.rc(max(xs) + pad, min(ys) - pad)
    r0, c0 = max(int(r0), 0), max(int(c0), 0)
    r1, c1 = min(int(r1), srf.S.shape[0] - 1), min(int(c1), srf.S.shape[1] - 1)
    S = srf.S[r0:r1, c0:c1]; G = srf.G[r0:r1, c0:c1]
    hb = S - G
    Tg = srf.Tg[r0:r1, c0:c1] if srf.Tg is not None else np.ones_like(S)
    # local maxima of the surface in a 150 m window
    mx = ndimage.maximum_filter(S, size=15)
    cand_mask = (S == mx) & (((hb >= args.min_roof) & (Tg > 0)) | (hb < 2))
    rr, cc = np.nonzero(cand_mask)
    order = np.argsort(-S[rr, cc])[: args.max_candidates]
    cands = []
    for k in order:
        r, c = rr[k] + r0, cc[k] + c0
        x = srf.T.c + (c + 0.5) * srf.T.a
        y = srf.T.f + (r + 0.5) * srf.T.e
        roof = hb[rr[k], cc[k]] >= 2
        z = float(srf.S[r, c]) + args.rep_mast
        cands.append(dict(name="R", x=x, y=y, z=z, base=float(srf.S[r, c]), roof=roof,
                          bh=float(hb[rr[k], cc[k]])))
    rep_budget_up = args.tx_dbm + args.gain_tx + args.rep_gain - args.losses   # handheld -> repeater
    rep_budget_dn = args.rep_tx_dbm + args.rep_gain + args.gain_rx - args.losses  # repeater -> handheld
    res = []
    for cd in cands:
        ms = []
        for p in pts:
            L = link(srf, cd, p, args.freq)
            m = min(rep_budget_up, rep_budget_dn) - L["loss"] - args.sens_dbm
            ms.append(m)
        ms = np.array(ms)
        res.append(dict(cd, margins=ms, served=int((ms >= args.serve_margin).sum()), worst=float(ms.min())))
    res.sort(key=lambda r: (-r["served"], -r["worst"]))
    # spread: skip candidates within 300 m of a better one
    picked = []
    for r in res:
        if all(math.hypot(r["x"] - q["x"], r["y"] - q["y"]) > 300 for q in picked):
            picked.append(r)
        if len(picked) >= top:
            break
    for r in picked:
        r["lon"], r["lat"] = srf.to_ll.transform(r["x"], r["y"])
    return picked


def write_repeaters_csv(picked, pts, path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["rank", "lat", "lon", "type", "building_height_m", "site_elev_m_asl",
                    "points_served", "worst_margin_dB"] + [f"margin_{p['name']}" for p in pts])
        for i, r in enumerate(picked, 1):
            w.writerow([i, f"{r['lat']:.6f}", f"{r['lon']:.6f}", site_type(r),
                        f"{r['bh']:.0f}", f"{r['base']:.0f}", f"{r['served']}/{len(pts)}",
                        f"{r['worst']:.1f}"] + [f"{m:.1f}" for m in r["margins"]])


def site_type(r):
    if not r["roof"]:
        return "ground mast"
    if r["bh"] >= 190:
        return "check: probably a tower/mast, not a roof"
    return "rooftop"


COLORS = {"good": "#1a9850", "marginal": "#fee08b", "weak": "#f46d43", "no link": "#a50026"}


def write_geo(srf, pts, rows, reps, args):
    byname = {p["name"]: p for p in pts}
    feats = []
    for p in pts:
        feats.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [p["lon"], p["lat"]]},
                      "properties": {"kind": "point", "name": p["name"], "antenna_m": p["h"],
                                     "on_roof": p["roof"], "elev_asl": round(p["base"], 1)}})
    for r in rows:
        a, b = byname[r["a"]], byname[r["b"]]
        feats.append({"type": "Feature", "geometry": {"type": "LineString",
                      "coordinates": [[a["lon"], a["lat"]], [b["lon"], b["lat"]]]},
                      "properties": {"kind": "link", "from": r["a"], "to": r["b"], "status": r["status"],
                                     "margin_dB": round(r["margin"], 1), "los": r["los"],
                                     "fresnel": round(r["fresnel"], 2), "dist_km": round(r["dist_km"], 2),
                                     "color": COLORS[r["status"]]}})
    for i, r in enumerate(reps, 1):
        feats.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [r["lon"], r["lat"]]},
                      "properties": {"kind": "repeater", "rank": i, "served": f"{r['served']}/{len(pts)}",
                                     "worst_margin_dB": round(r["worst"], 1),
                                     "type": site_type(r)}})
    gj = {"type": "FeatureCollection", "features": feats}
    json.dump(gj, open(os.path.join(args.out, "radio_links.geojson"), "w"), ensure_ascii=False)
    write_kml(gj, os.path.join(args.out, "radio_links.kml"))
    write_html(gj, os.path.join(args.out, "radio_map.html"), args)


def write_kml(gj, path):
    def kcol(h):  # #rrggbb -> aabbggrr
        return "ff" + h[5:7] + h[3:5] + h[1:3]
    out = ['<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document>']
    for f in gj["features"]:
        p = f["properties"]; g = f["geometry"]
        if g["type"] == "Point":
            nm = p.get("name") or f"Repeater #{p['rank']} ({p['served']})"
            out.append(f"<Placemark><name>{nm}</name><Point><coordinates>{g['coordinates'][0]},{g['coordinates'][1]},0</coordinates></Point></Placemark>")
        else:
            c = " ".join(f"{x},{y},0" for x, y in g["coordinates"])
            out.append(f"<Placemark><name>{p['from']} - {p['to']}: {p['status']} ({p['margin_dB']} dB)</name>"
                       f"<Style><LineStyle><color>{kcol(p['color'])}</color><width>3</width></LineStyle></Style>"
                       f"<LineString><coordinates>{c}</coordinates></LineString></Placemark>")
    out.append("</Document></kml>")
    open(path, "w", encoding="utf-8").write("".join(out))


def write_html(gj, path, args):
    html = """<!DOCTYPE html><html><head><meta charset="utf-8"><title>Radio links</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>html,body,#m{height:100%;margin:0;font:13px system-ui}
.lg{background:#fff;padding:8px 10px;border-radius:6px;box-shadow:0 1px 4px #0004;line-height:1.6}
.lg i{display:inline-block;width:18px;height:4px;margin-right:6px;vertical-align:middle}</style></head>
<body><div id="m"></div><script>
const gj=__GJ__;const info=__INFO__;
const m=L.map('m');L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',{maxZoom:19,attribution:'&copy; OpenStreetMap contributors'}).addTo(m);
const layer=L.geoJSON(gj,{
 style:f=>({color:f.properties.color,weight:f.properties.status==='no link'?2:4,dashArray:f.properties.los?null:'6 6',opacity:.9}),
 pointToLayer:(f,ll)=>f.properties.kind==='repeater'
   ?L.circleMarker(ll,{radius:8,color:'#5e3c99',fillColor:'#b2abd2',fillOpacity:.9,weight:2}).bindTooltip('R'+f.properties.rank,{permanent:f.properties.rank<=3,direction:'top'})
   :L.circleMarker(ll,{radius:7,color:'#000',fillColor:'#2b83ba',fillOpacity:1,weight:1.5}).bindTooltip(f.properties.name,{permanent:true,direction:'right'}),
 onEachFeature:(f,l)=>{const p=f.properties;l.bindPopup(p.kind==='link'
   ?`<b>${p.from} &harr; ${p.to}</b><br>${p.dist_km} km, ${p.status}<br>margin ${p.margin_dB} dB<br>line of sight: ${p.los?'yes':'no'}<br>Fresnel clearance: ${p.fresnel}`
   :p.kind==='repeater'?`<b>Repeater candidate #${p.rank}</b><br>${p.type}<br>serves ${p.served}<br>worst margin ${p.worst_margin_dB} dB`
   :`<b>${p.name}</b><br>antenna ${p.antenna_m} m ${p.on_roof?'above roof':'above ground'}<br>ground/roof ${p.elev_asl} m ASL`)}}).addTo(m);
m.fitBounds(layer.getBounds().pad(0.15));
const lg=L.control({position:'bottomleft'});lg.onAdd=()=>{const d=L.DomUtil.create('div','lg');
d.innerHTML=`<b>${info}</b><br><i style="background:#1a9850"></i>good (&ge;20 dB)<br><i style="background:#fee08b"></i>marginal (10&ndash;20)<br><i style="background:#f46d43"></i>weak (0&ndash;10)<br><i style="background:#a50026"></i>no link<br>dashed = no clear line of sight`;return d};lg.addTo(m);
</script></body></html>"""
    info = f"{args.freq} MHz, {args.tx_dbm} dBm, sens. {args.sens_dbm} dBm"
    open(path, "w", encoding="utf-8").write(html.replace("__GJ__", json.dumps(gj, ensure_ascii=False)).replace("__INFO__", json.dumps(info)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("points")
    md = "out" if os.path.exists("out/kyiv_surface.tif") else "model"
    ap.add_argument("--surface", default=f"{md}/kyiv_surface.tif")
    ap.add_argument("--ground", default=f"{md}/kyiv_ground.tif")
    ap.add_argument("--tagged", default=f"{md}/kyiv_height_is_tagged.tif")
    ap.add_argument("--freq", type=float, default=446.0, help="MHz")
    ap.add_argument("--tx-dbm", type=float, default=37.0, help="handheld TX power (37 dBm = 5 W)")
    ap.add_argument("--gain-tx", type=float, default=0.0, help="dBi")
    ap.add_argument("--gain-rx", type=float, default=0.0, help="dBi")
    ap.add_argument("--losses", type=float, default=3.0, help="cable/body losses, dB")
    ap.add_argument("--sens-dbm", type=float, default=-118.0, help="receiver sensitivity")
    ap.add_argument("--repeaters", action="store_true", help="search for repeater sites")
    ap.add_argument("--rep-tx-dbm", type=float, default=44.0, help="repeater TX power (44 dBm = 25 W)")
    ap.add_argument("--rep-gain", type=float, default=3.0, help="repeater antenna gain dBi")
    ap.add_argument("--rep-mast", type=float, default=3.0, help="repeater antenna height above roof/ground")
    ap.add_argument("--min-roof", type=float, default=25.0, help="min building height for rooftop candidates")
    ap.add_argument("--serve-margin", type=float, default=10.0, help="dB margin counted as 'served'")
    ap.add_argument("--search-pad", type=float, default=2000.0, help="search area padding, m")
    ap.add_argument("--max-candidates", type=int, default=600)
    ap.add_argument("--out", default="results")
    run(ap.parse_args())
