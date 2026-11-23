"""Export the surface model for the browser version (docs/data/), cut into tiles.

The model area is split into square tiles (default 10 km = 1000 x 1000 cells of 10 m).
The browser downloads only the tiles around the user's points.

Each tile is one gzip file `t_<col>_<row>.bin.gz` holding, little-endian, in this order:
  surface : int16 decimetres, TILE x TILE cells, row-wise delta coded
  ground  : int16 decimetres, (TILE/GF + 2)^2 cells at GF*10 m, centres at (k - 0.5) * GF cells
            from the tile corner (one cell of overlap all round, so bilinear sampling never
            needs a neighbouring tile), row-wise delta coded
  tagged  : 1 bit per surface cell (1 = building height comes from an OSM tag), bit order little
Tiles with nothing but nodata are skipped. meta.json describes the layout.
"""
import gzip, json, os, sys
import numpy as np
import rasterio
from radio_los import read_m

TILE = 1000   # cells per tile side (10 km at 10 m)
GF = 2        # ground grid = 20 m


def delta(a):
    d = a.astype(np.int32)
    d[:, 1:] = d[:, 1:] - d[:, :-1]
    return d.astype("<i2")


def main(src, dst="docs/data"):
    os.makedirs(dst, exist_ok=True)
    for f in os.listdir(dst):
        if f.startswith("t_") or f in ("surface.i16.gz", "ground30.i16.gz", "tagged.bits.gz"):
            os.remove(os.path.join(dst, f))
    ds = rasterio.open(f"{src}/kyiv_surface.tif")
    S = read_m(ds)
    G = read_m(rasterio.open(f"{src}/kyiv_ground.tif"))
    Tg = rasterio.open(f"{src}/kyiv_height_is_tagged.tif").read(1) > 0
    H, W = S.shape
    assert H % TILE == 0 and W % TILE == 0, "model size must be a multiple of the tile size"
    ny, nx = H // TILE, W // TILE
    # ground at GF x coarser grid, padded by one coarse cell all round
    Gc = G.reshape(H // GF, GF, W // GF, GF).mean(axis=(1, 3))
    Gp = np.pad(Gc, 1, mode="edge")
    g = TILE // GF
    tiles, total = {}, 0
    for ty in range(ny):
        for tx in range(nx):
            s = S[ty * TILE:(ty + 1) * TILE, tx * TILE:(tx + 1) * TILE]
            gr = Gp[ty * g:(ty + 1) * g + 2, tx * g:(tx + 1) * g + 2]
            tg = Tg[ty * TILE:(ty + 1) * TILE, tx * TILE:(tx + 1) * TILE]
            raw = delta(np.round(s * 10)).tobytes() + delta(np.round(gr * 10)).tobytes() + \
                np.packbits(tg.ravel(), bitorder="little").tobytes()
            name = f"t_{tx}_{ty}.bin.gz"
            with gzip.open(os.path.join(dst, name), "wb", compresslevel=9) as f:
                f.write(raw)
            n = os.path.getsize(os.path.join(dst, name))
            tiles[f"{tx}_{ty}"] = n
            total += n
    T = ds.transform
    meta = dict(version=2, crs="EPSG:32636", utm_zone=36, res=T.a, x0=T.c, y0=T.f, width=W, height=H,
                tile=TILE, ground_factor=GF, nx=nx, ny=ny, scale=0.1,
                layout="per tile gzip: surface int16 dm delta | ground int16 dm delta (TILE/GF+2)^2 | tagged bits",
                tiles=tiles, total_bytes=total)
    json.dump(meta, open(os.path.join(dst, "meta.json"), "w"))
    print(f"{nx}x{ny} tiles, {total / 1e6:.1f} MB total, largest {max(tiles.values()) / 1e6:.2f} MB")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ("out" if os.path.exists("out/kyiv_surface.tif") else "model"))
