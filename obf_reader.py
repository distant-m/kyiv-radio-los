"""Minimal reader for OsmAnd .obf map files (map index section only).

Follows the OBF.proto spec from osmandapp/OsmAnd-resources and the logic of
OsmAnd's BinaryMapIndexReader.java. Extracts map objects (with their tags and
coordinates) from the most detailed zoom level.
"""
import mmap
import math
import sys


class Buf:
    def __init__(self, mm, pos=0):
        self.mm = mm
        self.pos = pos

    def varint(self):
        mm, p = self.mm, self.pos
        result = 0
        shift = 0
        while True:
            b = mm[p]
            p += 1
            result |= (b & 0x7F) << shift
            if b < 0x80:
                break
            shift += 7
        self.pos = p
        return result

    def sint(self):
        v = self.varint()
        return (v >> 1) ^ -(v & 1)

    def int32be(self):
        p = self.pos
        v = int.from_bytes(self.mm[p:p + 4], "big")
        self.pos = p + 4
        return v

    def tag(self):
        t = self.varint()
        return t >> 3, t & 7

    def skip(self, wt):
        if wt == 0:
            self.varint()
        elif wt == 1:
            self.pos += 8
        elif wt == 2:
            n = self.varint()
            self.pos += n
        elif wt == 5:
            self.pos += 4
        elif wt == 6:  # OsmAnd: fixed32 big-endian length delimited
            n = self.int32be()
            self.pos += n
        else:
            raise ValueError(f"unknown wire type {wt} at {self.pos}")

    def string(self):
        n = self.varint()
        s = self.mm[self.pos:self.pos + n].decode("utf-8", "replace")
        self.pos += n
        return s


def x31_to_lon(x):
    return x / 2147483648.0 * 360.0 - 180.0


def y31_to_lat(y):
    n = math.pi - 2.0 * math.pi * y / 2147483648.0
    return math.degrees(math.atan(math.sinh(n)))


def lon_to_x31(lon):
    return int((lon + 180.0) / 360.0 * 2147483648.0)


def lat_to_y31(lat):
    r = math.radians(lat)
    return int((1.0 - math.log(math.tan(r) + 1.0 / math.cos(r)) / math.pi) / 2.0 * 2147483648.0)


class OBF:
    def __init__(self, path):
        self.f = open(path, "rb")
        self.mm = mmap.mmap(self.f.fileno(), 0, access=mmap.ACCESS_READ)
        self.map_indexes = []  # (start, end)
        b = Buf(self.mm)
        end = len(self.mm)
        while b.pos < end:
            field, wt = b.tag()
            if field == 0:
                break
            if field == 6 and wt == 6:
                ln = b.int32be()
                self.map_indexes.append((b.pos, b.pos + ln))
                b.pos += ln
            elif field == 32:
                b.varint()
                break
            else:
                b.skip(wt)

    def read_map_index(self, idx=0):
        start, end = self.map_indexes[idx]
        b = Buf(self.mm, start)
        rules = {}
        levels = []
        name = None
        next_id = 1
        while b.pos < end:
            field, wt = b.tag()
            if field == 2:
                name = b.string()
            elif field == 4:
                ln = b.varint()
                rend = b.pos + ln
                tag = value = None
                rid = None
                rtype = 0
                while b.pos < rend:
                    f2, w2 = b.tag()
                    if f2 == 3:
                        tag = b.string()
                    elif f2 == 5:
                        value = b.string()
                    elif f2 == 7:
                        rid = b.varint()
                    elif f2 == 10:
                        rtype = b.varint()
                    else:
                        b.skip(w2)
                if rid is None:
                    rid = next_id
                next_id = rid + 1
                rules[rid] = (tag, value, rtype)
            elif field == 5 and wt == 6:
                ln = b.int32be()
                lstart = b.pos
                lv = {"start": lstart, "end": lstart + ln}
                # read header fields until boxes
                lb = Buf(self.mm, lstart)
                while lb.pos < lstart + ln:
                    f2, w2 = lb.tag()
                    if f2 in (1, 2, 3, 4, 5, 6):
                        lv[{1: "maxZoom", 2: "minZoom", 3: "left", 4: "right",
                            5: "top", 6: "bottom"}[f2]] = lb.varint()
                    else:
                        break
                levels.append(lv)
                b.pos += ln
            else:
                b.skip(wt)
        self.name, self.rules, self.levels = name, rules, levels
        return name, rules, levels

    def iter_leaf_boxes(self, level, bbox31=None):
        """Yield (left, top, mapDataBlockOffset) for boxes carrying data."""
        out = []
        b = Buf(self.mm, level["start"])
        end = level["end"]
        root = dict(left=level["left"], right=level["right"], top=level["top"], bottom=level["bottom"])
        while b.pos < end:
            field, wt = b.tag()
            if field == 7 and wt == 6:
                ln = b.int32be()
                fp = b.pos
                self._read_box(fp, fp + ln, root, bbox31, out)
                b.pos = fp + ln
            else:
                b.skip(wt)
        return out

    def _read_box(self, start, end, parent, bbox31, out):
        b = Buf(self.mm, start)
        cur = {}
        while b.pos < end:
            field, wt = b.tag()
            if field == 1:
                cur["left"] = b.sint() + parent["left"]
            elif field == 2:
                cur["right"] = b.sint() + parent["right"]
            elif field == 3:
                cur["top"] = b.sint() + parent["top"]
            elif field == 4:
                cur["bottom"] = b.sint() + parent["bottom"]
                if bbox31:
                    l, t, r, bo = bbox31
                    if cur["right"] < l or cur["left"] > r or cur["bottom"] < t or cur["top"] > bo:
                        return
            elif field == 5:
                out.append((cur["left"], cur["top"], b.int32be() + start))
            elif field == 6:
                b.varint()
            elif field == 7:
                ln = b.int32be()
                fp = b.pos
                self._read_box(fp, fp + ln, cur, bbox31, out)
                b.pos = fp + ln
            else:
                b.skip(wt)

    def read_block(self, left, top, offset, want=None):
        """Read a MapDataBlock. want(types, addtypes) -> bool filters objects."""
        b = Buf(self.mm, offset)
        ln = b.varint()
        end = b.pos + ln
        objs = []
        strings = []
        while b.pos < end:
            field, wt = b.tag()
            if field == 10:
                b.varint()
            elif field == 12:
                l2 = b.varint()
                o = self._read_obj(b.pos, b.pos + l2, left, top, want)
                b.pos += l2
                if o is not None:
                    objs.append(o)
            elif field == 15:
                l2 = b.varint()
                send = b.pos + l2
                while b.pos < send:
                    f2, w2 = b.tag()
                    if f2 == 1:
                        strings.append(b.string())
                    else:
                        b.skip(w2)
            else:
                b.skip(wt)
        for o in objs:
            names = o.pop("_names", None)
            if names:
                o["strings"] = {k: strings[v] if v < len(strings) else None for k, v in names}
        return objs

    def _coords(self, b, end, left, top):
        mask = ~((1 << 5) - 1)
        px, py = left & mask, top & mask
        pts = []
        while b.pos < end:
            x = (b.sint() << 5) + px
            y = (b.sint() << 5) + py
            pts.append((x, y))
            px, py = x, y
        return pts

    def _read_obj(self, start, end, left, top, want):
        b = Buf(self.mm, start)
        o = {}
        coord_span = None
        while b.pos < end:
            field, wt = b.tag()
            if field in (1, 2):
                ln = b.varint()
                coord_span = (b.pos, b.pos + ln)
                o["area"] = field == 2
                b.pos += ln
            elif field == 4:
                ln = b.varint()
                o.setdefault("_inner_spans", []).append((b.pos, b.pos + ln))
                b.pos += ln
            elif field in (6, 7):
                ln = b.varint()
                e2 = b.pos + ln
                lst = []
                while b.pos < e2:
                    lst.append(b.varint())
                o["types" if field == 7 else "addtypes"] = lst
            elif field == 10:
                ln = b.varint()
                e2 = b.pos + ln
                lst = []
                while b.pos < e2:
                    lst.append((b.varint(), b.varint()))
                o["_names"] = lst
            elif field == 12:
                o["id"] = b.sint()
            else:
                b.skip(wt)
        if want is not None and not want(o.get("types", []), o.get("addtypes", [])):
            return None
        if coord_span:
            o["coords"] = self._coords(Buf(self.mm, coord_span[0]), coord_span[1], left, top)
        inner = o.pop("_inner_spans", [])
        if inner:
            o["inner"] = [self._coords(Buf(self.mm, s), e, left, top) for s, e in inner]
        return o


if __name__ == "__main__":
    obf = OBF(sys.argv[1])
    name, rules, levels = obf.read_map_index()
    print("map index:", name, "rules:", len(rules))
    for lv in levels:
        print("level", lv["minZoom"], "-", lv["maxZoom"], "size MB", round((lv["end"] - lv["start"]) / 1e6, 1))
    kw = ("build", "height", "level", "roof", "min_")
    for rid, (t, v, ty) in sorted(rules.items()):
        if any(k in (t or "") for k in kw):
            print(rid, t, v, ty)
