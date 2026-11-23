"""Extract building footprints + heights from an OsmAnd .obf into GeoJSON-like records."""
import json, sys, time
from obf_reader import OBF, x31_to_lon, y31_to_lat

def main(path, out):
    obf = OBF(path)
    name, rules, levels = obf.read_map_index()
    lv = max(levels, key=lambda l: l["maxZoom"])
    bld_ids = {rid for rid, (t, v, ty) in rules.items() if t in ("building", "building:part")}
    text_ids = {rid: t for rid, (t, v, ty) in rules.items()
                if t in ("building:levels", "height", "min_height", "building:min_level", "roof:levels", "roof:height", "ele")}
    want = lambda types, add: any(t in bld_ids for t in types)
    t0 = time.time()
    boxes = obf.iter_leaf_boxes(lv)
    seen = set(); feats = []
    for left, top, off in boxes:
        for o in obf.read_block(left, top, off, want):
            if o["id"] in seen or "coords" not in o:
                continue
            seen.add(o["id"])
            tags = {}
            for t in o.get("types", []) + o.get("addtypes", []):
                tg, v, _ = rules[t]
                if v is not None:
                    tags[tg] = v
            for k, v in (o.get("strings") or {}).items():
                if k in text_ids:
                    tags[text_ids[k]] = v
            ring = [(round(x31_to_lon(x), 7), round(y31_to_lat(y), 7)) for x, y in o["coords"]]
            feats.append({"id": o["id"], "tags": tags, "ring": ring})
    print(f"boxes {len(boxes)} buildings {len(feats)} in {time.time()-t0:.0f}s")
    json.dump(feats, open(out, "w"))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
