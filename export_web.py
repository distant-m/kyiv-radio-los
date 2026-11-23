"""Export the surface model for the browser version (docs/).

surface: int16 decimetres, 10 m grid, row-wise delta coded
ground:  int16 decimetres, 30 m grid (terrain is 30 m data anyway; sampled bilinearly)
tagged:  1 bit per 10 m cell (1 = building height comes from an OSM tag)
All gzip-compressed little-endian, described in meta.json.
"""
import gzip, json, os, sys
import numpy as np
import rasterio
from radio_los import read_m

src = sys.argv[1] if len(sys.argv) > 1 else ("out" if os.path.exists("out/kyiv_surface.tif") else "model")
dst = "docs/data"
os.makedirs(dst, exist_ok=True)
ds = rasterio.open(f"{src}/kyiv_surface.tif")
S = read_m(ds); G = read_m(rasterio.open(f"{src}/kyiv_ground.tif"))
Tg = rasterio.open(f"{src}/kyiv_height_is_tagged.tif").read(1) > 0
H, W = S.shape
f = 3
H3, W3 = (H + f - 1) // f, (W + f - 1) // f
Gp = np.pad(G, ((0, H3 * f - H), (0, W3 * f - W)), mode="edge")
G3 = Gp.reshape(H3, f, W3, f).mean(axis=(1, 3))

def delta(a):
    """Row-wise difference (first column kept) - compresses ~2x better; JS does a running sum."""
    d = a.astype(np.int32)
    d[:, 1:] = d[:, 1:] - d[:, :-1]
    return d.astype("<i2")


def save(name, arr):
    with gzip.open(f"{dst}/{name}", "wb", compresslevel=9) as g:
        g.write(arr.tobytes())
    return os.path.getsize(f"{dst}/{name}")

sizes = {
    "surface": save("surface.i16.gz", delta(np.round(S * 10))),
    "ground": save("ground30.i16.gz", delta(np.round(G3 * 10))),
    "tagged": save("tagged.bits.gz", np.packbits(Tg.ravel(), bitorder="little")),
}
T = ds.transform
meta = dict(crs="EPSG:32636", utm_zone=36, width=W, height=H, x0=T.c, y0=T.f, res=T.a,
            encoding="int16 decimetres, row-wise delta, gzip", ground_factor=f, ground_width=W3, ground_height=H3, scale=0.1,
            files=dict(surface="surface.i16.gz", ground="ground30.i16.gz", tagged="tagged.bits.gz"),
            bytes=sizes)
json.dump(meta, open(f"{dst}/meta.json", "w"), indent=1)
print(meta)
