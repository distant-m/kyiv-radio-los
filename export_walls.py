"""Export building outlines (outer walls) per web tile, for the reflection model.

Source: buildings.json from extract_buildings.py (OsmAnd .obf). Only whole buildings
(everything except `building:part` outlines) are used; outlines are simplified (0.5 m),
oriented counter-clockwise (outside is to the right of each edge), and buildings
under 40 m2 are skipped. Each building goes to the tile holding its centroid.

File docs/data/w_<col>_<row>.bin.gz = int16 little-endian stream, per building:
  n, x0, y0, dx1, dy1, ... dx(n-1), dy(n-1)
in 0.5 m units; x0/y0 relative to the tile's south-west corner (can be slightly outside).
Wall heights are not stored: the browser reads them from the surface model.
"""
import gzip, json, os, sys
import numpy as np
from pyproj import Transformer
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
from shapely import transform as stf

def main(buildings_json, dst="docs/data"):
    meta = json.load(open(os.path.join(dst, "meta.json")))
    T, res, x0, y0top = meta["tile"], meta["res"], meta["x0"], meta["y0"]
    tm = T * res
    nx, ny = meta["nx"], meta["ny"]
    to_utm = Transformer.from_crs("EPSG:4326", meta["crs"], always_xy=True)
    per = {}
    kept = skipped = 0
    for f in json.load(open(buildings_json)):
        if ("building:part" in f["tags"] and "building" not in f["tags"]) or len(f["ring"]) < 4:
            continue
        ring = np.array(f["ring"])
        xs, ys = to_utm.transform(ring[:, 0], ring[:, 1])
        if xs.max() < x0 or xs.min() > x0 + nx * tm or ys.max() < y0top - ny * tm or ys.min() > y0top:
            continue
        try:
            p = Polygon(np.column_stack([xs, ys]))
            if not p.is_valid:
                p = p.buffer(0)
            if p.geom_type == "MultiPolygon":
                p = max(p.geoms, key=lambda g: g.area)
            if p.is_empty or p.area < 40:
                skipped += 1
                continue
            p = orient(p.simplify(0.5, preserve_topology=True), 1.0)
        except Exception:
            continue
        cx, cy = p.centroid.x, p.centroid.y
        tx, ty = int((cx - x0) // tm), int((y0top - cy) // tm)
        if not (0 <= tx < nx and 0 <= ty < ny):
            continue
        ox, oy = x0 + tx * tm, y0top - (ty + 1) * tm   # tile SW corner
        pts = np.round((np.array(p.exterior.coords)[:-1] - [ox, oy]) * 2).astype(np.int32)
        if len(pts) < 3 or np.abs(pts).max() > 32000:
            continue
        d = np.vstack([pts[:1], np.diff(pts, axis=0)])
        if np.abs(d[1:]).max(initial=0) > 32000:
            continue
        per.setdefault(f"{tx}_{ty}", []).append(np.concatenate([[len(pts)], d.ravel()]))
        kept += 1
    total = 0
    walls = {}
    for k, arrs in per.items():
        raw = np.concatenate(arrs).astype("<i2").tobytes()
        name = os.path.join(dst, f"w_{k}.bin.gz")
        with gzip.open(name, "wb", compresslevel=9) as g:
            g.write(raw)
        walls[k] = os.path.getsize(name)
        total += walls[k]
    meta["walls"] = walls
    meta["walls_info"] = "per tile gzip int16 stream: n, x0, y0, dx, dy ... (0.5 m units, from tile SW corner), CCW outer rings"
    json.dump(meta, open(os.path.join(dst, "meta.json"), "w"))
    print(f"buildings kept {kept}, skipped small {skipped}; wall files {len(walls)}, {total / 1e6:.1f} MB")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/buildings.json")
