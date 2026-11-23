"""Local map app for kyiv-radio-los.

    python3 app.py            # then open http://localhost:8000

Uses only the Python standard library for the web server. Needs the surface
rasters built by build_surface.py (default: out/kyiv_*.tif).
"""
import argparse, json, math, os, sys, traceback, webbrowser
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from types import SimpleNamespace

import radio_los as R

HERE = os.path.dirname(os.path.abspath(__file__))

DEFAULTS = dict(freq=446.0, tx_dbm=37.0, sens_dbm=-118.0, gain_tx=0.0, gain_rx=0.0, losses=3.0,
                repeaters=False, rep_tx_dbm=44.0, rep_gain=3.0, rep_mast=3.0, min_roof=25.0,
                serve_margin=10.0, search_pad=2000.0, max_candidates=600)
LIMITS = dict(freq=(30, 6000), tx_dbm=(-10, 60), sens_dbm=(-150, -60), gain_tx=(-10, 30), gain_rx=(-10, 30),
              losses=(0, 40), rep_tx_dbm=(0, 60), rep_gain=(-10, 30), rep_mast=(0, 100), min_roof=(0, 300),
              serve_margin=(-20, 60), search_pad=(0, 10000), max_candidates=(10, 3000))


class Planner:
    def __init__(self, surface, ground, tagged):
        self.srf = R.Surface(surface, ground, tagged)
        b = self.srf.ds.bounds
        lons, lats = self.srf.to_ll.transform([b.left, b.right, b.left, b.right], [b.bottom, b.bottom, b.top, b.top])
        self.bounds = [[max(lats[:2]), max(lons[0], lons[2])], [min(lats[2:]), min(lons[1], lons[3])]]

    def params(self, p):
        out = dict(DEFAULTS)
        for k, v in (p or {}).items():
            if k == "repeaters":
                out[k] = bool(v)
            elif k in LIMITS and v not in (None, ""):
                lo, hi = LIMITS[k]
                v = float(v)
                if not (lo <= v <= hi) or math.isnan(v):
                    raise ValueError(f"{k} = {v} is outside {lo}..{hi}")
                out[k] = int(v) if k == "max_candidates" else v
        return SimpleNamespace(**out)

    def place_points(self, raw):
        pts, names = [], set()
        for i, p in enumerate(raw):
            name = str(p.get("name") or f"P{i + 1}").strip()[:40]
            if name in names:
                name = f"{name} ({i + 1})"
            names.add(name)
            q = dict(name=name, lat=float(p["lat"]), lon=float(p["lon"]),
                     h=float(p.get("h", 1.5)), roof=bool(p.get("roof", False)))
            if not (0 <= q["h"] <= 500):
                raise ValueError(f"{name}: antenna height {q['h']} m is not plausible")
            x, y = self.srf.to_xy.transform(q["lon"], q["lat"])
            if not self.srf.inside(x, y):
                raise ValueError(f"{name} is outside the map area of the surface model")
            pts.append(R.place(self.srf, q))
        return pts

    def elevation(self, lat, lon):
        x, y = self.srf.to_xy.transform(lon, lat)
        if not self.srf.inside(x, y):
            return dict(inside=False)
        g = float(self.srf.sample(self.srf.G, x, y)); s = float(self.srf.sample(self.srf.S, x, y))
        tg = bool(self.srf.Tg is not None and self.srf.sample(self.srf.Tg, x, y) > 0)
        return dict(inside=True, ground=round(g, 1), surface=round(s, 1), building=round(s - g, 1) if s - g > 1 else 0,
                    height_tagged=tg)

    def analyze(self, body):
        args = self.params(body.get("params"))
        pts = self.place_points(body.get("points", []))
        if len(pts) < 1:
            raise ValueError("Add at least one point")
        rows, budget = R.compute_links(self.srf, pts, args) if len(pts) >= 2 else ([], 0)
        reps = R.repeater_search(self.srf, pts, args, budget) if args.repeaters else []
        return dict(
            points=[dict(name=p["name"], lat=p["lat"], lon=p["lon"], h=p["h"], roof=p["roof"],
                         base=round(p["base"], 1), antenna_asl=round(p["z"], 1)) for p in pts],
            links=[dict(a=r["a"], b=r["b"], dist_km=round(r["dist_km"], 2), los=bool(r["los"]),
                        fresnel=round(r["fresnel"], 2), fspl=round(r["fspl"], 1), diff=round(r["diff"], 1),
                        model=r["model"], rx_dbm=round(r["rx_dbm"], 1), margin=round(r["margin"], 1),
                        status=r["status"], obstacle="building" if r["obstacle_is_building"] else "terrain",
                        obstacle_height="tagged" if r["obstacle_height_known"] else "estimated")
                   for r in rows],
            repeaters=[dict(rank=i, lat=round(r["lat"], 6), lon=round(r["lon"], 6), type=R.site_type(r),
                            building_m=round(r["bh"]), site_asl=round(r["base"]), served=r["served"],
                            total=len(pts), worst=round(r["worst"], 1),
                            margins={p["name"]: round(float(m), 1) for p, m in zip(pts, r["margins"])})
                       for i, r in enumerate(reps, 1)],
            params=vars(args))

    def profile(self, body):
        args = self.params(body.get("params"))
        a, b = self.place_points([body["a"], body["b"]])
        pr = R.profile(self.srf, a, b, args.freq)
        L = R.link(self.srf, a, b, args.freq)
        pr.update(names=[a["name"], b["name"]], los=bool(L["los"]), fresnel=round(L["fresnel"], 2),
                  dist_km=round(L["dist_km"], 2), freq=args.freq)
        return pr


class Handler(BaseHTTPRequestHandler):
    planner = None

    def log_message(self, fmt, *a):
        pass

    def send(self, code, body, ctype="application/json"):
        data = body if isinstance(body, bytes) else json.dumps(body, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype + ("; charset=utf-8" if "json" in ctype or "html" in ctype else ""))
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path in ("/", "/index.html"):
            with open(os.path.join(HERE, "static", "index.html"), "rb") as f:
                return self.send(200, f.read(), "text/html")
        if self.path == "/api/info":
            return self.send(200, dict(bounds=self.planner.bounds, defaults=DEFAULTS))
        self.send(404, {"error": "not found"})

    def do_POST(self):
        try:
            n = int(self.headers.get("Content-Length", 0))
            if n > 2_000_000:
                return self.send(413, {"error": "request too large"})
            body = json.loads(self.rfile.read(n) or b"{}")
            if self.path == "/api/analyze":
                return self.send(200, self.planner.analyze(body))
            if self.path == "/api/profile":
                return self.send(200, self.planner.profile(body))
            if self.path == "/api/elevation":
                return self.send(200, self.planner.elevation(float(body["lat"]), float(body["lon"])))
            self.send(404, {"error": "not found"})
        except (ValueError, KeyError, TypeError) as e:
            self.send(400, {"error": str(e)})
        except Exception as e:
            traceback.print_exc()
            self.send(500, {"error": f"internal error: {e}"})


def main():
    ap = argparse.ArgumentParser()
    # A model you built yourself (out/) wins over the ready-made one shipped in the repo (model/).
    md = os.path.join(HERE, "out") if os.path.exists(os.path.join(HERE, "out", "kyiv_surface.tif")) else os.path.join(HERE, "model")
    ap.add_argument("--surface", default=os.path.join(md, "kyiv_surface.tif"))
    ap.add_argument("--ground", default=os.path.join(md, "kyiv_ground.tif"))
    ap.add_argument("--tagged", default=os.path.join(md, "kyiv_height_is_tagged.tif"))
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args()
    for f in (a.surface, a.ground):
        if not os.path.exists(f):
            sys.exit(f"Missing {f}. Build it first:\n  python3 build_surface.py data/N50E029.tif data/N50E030.tif data/buildings.json out/kyiv")
    print(f"Loading surface model from {os.path.dirname(a.surface)} ...")
    Handler.planner = Planner(a.surface, a.ground, a.tagged)
    srv = ThreadingHTTPServer(("127.0.0.1", a.port), Handler)
    url = f"http://localhost:{a.port}"
    print(f"Radio planner running at {url}  (Ctrl+C to stop)")
    if not a.no_browser:
        webbrowser.open(url)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
