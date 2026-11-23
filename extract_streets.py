"""Extract the street network (OSM highway ways) from an OsmAnd .obf map.

Output JSON: list of {"hw": highway type, "pts": [[x31, y31], ...]} using OsmAnd's raw
31-bit tile coordinates, so a junction shared by two streets has identical coordinates.
Point features (lamps, crossings...), tunnels, tracks, bridleways, raceways and proposed
roads are skipped.
"""
import json, sys, time
from obf_reader import OBF

KEEP = {"footway", "service", "residential", "path", "unclassified", "tertiary", "secondary", "steps",
        "primary", "living_street", "trunk", "trunk_link", "cycleway", "primary_link", "secondary_link",
        "pedestrian", "tertiary_link", "motorway", "motorway_link", "road"}

def main(path, out):
    obf = OBF(path)
    name, rules, levels = obf.read_map_index()
    lv = max(levels, key=lambda l: l["maxZoom"])
    hw = {rid: v for rid, (t, v, ty) in rules.items() if t == "highway" and v in KEEP}
    tunnel = {rid for rid, (t, v, ty) in rules.items() if t == "tunnel" and v == "yes"}
    want = lambda types, add: any(t in hw for t in types)
    t0 = time.time(); seen = set(); ways = []
    for left, top, off in obf.iter_leaf_boxes(lv):
        for o in obf.read_block(left, top, off, want):
            if o["id"] in seen or o.get("area") or "coords" not in o or len(o["coords"]) < 2:
                continue
            seen.add(o["id"])
            if any(t in tunnel for t in o.get("types", []) + o.get("addtypes", [])):
                continue
            h = next(hw[t] for t in o["types"] if t in hw)
            ways.append({"hw": h, "pts": [list(p) for p in o["coords"]]})
    print(f"{len(ways)} street ways in {time.time() - t0:.0f}s")
    json.dump(ways, open(out, "w"))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
