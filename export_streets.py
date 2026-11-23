"""Export the street network per web tile, for the street-corridor model.

Input: streets.json from extract_streets.py (OsmAnd raw coordinates, junctions shared).
Ways are split at junctions (vertices used by more than one way, or twice by one way).
Loose ends within 6 m of another street are joined to it (OSM driveways and paths often stop
just short of the road), splitting that street at the joint. Each piece is simplified (Douglas-Peucker, 6 m: a street that stays inside its own width
counts as straight), and pieces longer than 4 km are split. A piece goes to the tile of
its first vertex.

File docs/data/s_<col>_<row>.bin.gz = int16 little-endian stream, per piece:
  n, x0, y0, dx1, dy1, ... (0.5 m units, x0/y0 from the tile's south-west corner)
Pieces meet at identical end coordinates, which is how the browser rebuilds the graph.
"""
import gzip, json, os, sys, math
from collections import Counter
import numpy as np
from pyproj import Transformer
from shapely.geometry import LineString, Point
from shapely import STRtree

def x31_lon(x): return x / 2147483648.0 * 360.0 - 180.0
def y31_lat(y):
    n = math.pi - 2.0 * math.pi * y / 2147483648.0
    return math.degrees(math.atan(math.sinh(n)))

def main(streets_json, dst="docs/data"):
    meta = json.load(open(os.path.join(dst, "meta.json")))
    T, res, X0, Y0 = meta["tile"], meta["res"], meta["x0"], meta["y0"]
    tm = T * res; nx, ny = meta["nx"], meta["ny"]
    xmin, xmax, ymin, ymax = X0 - 2000, X0 + nx * tm + 2000, Y0 - ny * tm - 2000, Y0 + 2000
    to_utm = Transformer.from_crs("EPSG:4326", meta["crs"], always_xy=True)
    ways = json.load(open(streets_json))
    use = Counter()
    for w in ways:
        for p in w["pts"]:
            use[(p[0], p[1])] += 1
    pieces = []
    for w in ways:
        pts = [tuple(p) for p in w["pts"]]
        cut = [0] + [i for i in range(1, len(pts) - 1) if use[pts[i]] > 1] + [len(pts) - 1]
        for a, b in zip(cut[:-1], cut[1:]):
            if b > a:
                pieces.append(pts[a:b + 1])
    # 1. pieces in UTM, snapped to the 0.5 m grid
    lines = []
    for pc in pieces:
        lon = np.array([x31_lon(p[0]) for p in pc]); lat = np.array([y31_lat(p[1]) for p in pc])
        xs, ys = to_utm.transform(lon, lat)
        if xs.max() < xmin or xs.min() > xmax or ys.max() < ymin or ys.min() > ymax:
            continue
        xy = np.round(np.column_stack([xs, ys]) * 2) / 2
        if len(xy) >= 2 and np.hypot(*(xy[-1] - xy[0])) + np.abs(np.diff(xy, axis=0)).sum() > 0:
            lines.append(xy)
    # 2. repair: a loose end within SNAP m of another street joins it there
    #    (OSM driveways and paths often stop just short of the road or meet it without a shared node)
    SNAP = 6.0
    ends = Counter()
    for xy in lines:
        ends[tuple(xy[0])] += 1; ends[tuple(xy[-1])] += 1
    geoms = [LineString(xy) for xy in lines]
    tree = STRtree(geoms)
    splits = {}   # piece index -> list of (distance along, point)
    moved = 0
    for i, xy in enumerate(lines):
        for end in (0, -1):
            p = tuple(xy[end])
            if ends[p] > 1:
                continue
            pt = Point(p)
            best = None
            for j in tree.query(pt.buffer(SNAP)):
                if j == i:
                    continue
                d = geoms[j].distance(pt)
                if d <= SNAP and (best is None or d < best[0]):
                    best = (d, j)
            if best is None:
                continue
            j = best[1]
            along = geoms[j].project(pt)
            q = geoms[j].interpolate(along)
            q = (round(q.x * 2) / 2, round(q.y * 2) / 2)
            lines[i][end] = q
            splits.setdefault(j, []).append((along, q))
            moved += 1
    final = []
    for j, xy in enumerate(lines):
        if j not in splits:
            final.append(xy); continue
        cuts = sorted(set(splits[j]))
        g = LineString(xy); prev = 0.0; start = [tuple(xy[0])]
        coords = [tuple(c) for c in xy]
        # walk the polyline, inserting split points in order
        acc = [0.0]
        for k in range(1, len(coords)):
            acc.append(acc[-1] + float(np.hypot(coords[k][0] - coords[k - 1][0], coords[k][1] - coords[k - 1][1])))
        cur = [coords[0]]; ci = 0
        for k in range(1, len(coords)):
            while ci < len(cuts) and cuts[ci][0] <= acc[k]:
                q = cuts[ci][1]
                if q != cur[-1]:
                    cur.append(q)
                if len(cur) > 1:
                    final.append(np.array(cur))
                cur = [q]; ci += 1
            if coords[k] != cur[-1]:
                cur.append(coords[k])
        if len(cur) > 1:
            final.append(np.array(cur))
    print(f"loose ends joined to a nearby street: {moved}")
    out = {}
    kept = 0
    for xy in final:
        if len(xy) > 2:
            xy = np.array(LineString(xy).simplify(6.0, preserve_topology=False).coords)
        if len(xy) < 2 or (len(xy) == 2 and np.allclose(xy[0], xy[1])):
            continue
        # split long pieces (keeps int16 offsets small)
        chunks, cur, length = [], [xy[0]], 0.0
        for p in xy[1:]:
            length += float(np.hypot(*(p - cur[-1]))); cur.append(p)
            if length > 4000:
                chunks.append(cur); cur, length = [p], 0.0
        if len(cur) > 1:
            chunks.append(cur)
        for ch in chunks:
            ch = np.array(ch)
            tx = int((ch[0, 0] - X0) // tm); ty = int((Y0 - ch[0, 1]) // tm)
            tx = min(max(tx, 0), nx - 1); ty = min(max(ty, 0), ny - 1)
            ox, oy = X0 + tx * tm, Y0 - (ty + 1) * tm
            q = np.round((ch - [ox, oy]) * 2).astype(np.int64)
            d = np.vstack([q[:1], np.diff(q, axis=0)])
            if np.abs(d).max() > 32000:
                continue
            out.setdefault(f"{tx}_{ty}", []).append(np.concatenate([[len(q)], d.ravel()]))
            kept += 1
    files, total = {}, 0
    for k, arrs in out.items():
        name = os.path.join(dst, f"s_{k}.bin.gz")
        with gzip.open(name, "wb", compresslevel=9) as g:
            g.write(np.concatenate(arrs).astype("<i2").tobytes())
        files[k] = os.path.getsize(name); total += files[k]
    meta["streets"] = files
    meta["streets_info"] = "per tile gzip int16 stream: n, x0, y0, dx, dy ... (0.5 m units, from tile SW corner), open polylines meeting at identical ends"
    json.dump(meta, open(os.path.join(dst, "meta.json"), "w"))
    print(f"street pieces {kept} in {len(files)} tiles, {total / 1e6:.1f} MB")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/streets.json")
