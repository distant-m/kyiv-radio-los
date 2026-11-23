# kyiv-radio-los

Radio visibility and repeater placement for many points at once, using **terrain + buildings**.

For a set of radio points (handhelds, base stations) it answers:

1. **Which pairs can talk directly?** It reports line of sight, Fresnel-zone clearance, estimated path loss and link margin for every pair.
2. **Where should a repeater go?** It ranks rooftops and high ground by how many of your points each one reaches, and with what margin.

The results come as CSV tables, a KML file (Google Earth, OsmAnd, phone map apps), GeoJSON (QGIS), and an HTML map you can open in any browser.

## Quick start

The repo already contains a ready-made Kyiv model (terrain + buildings, `model/`, 16 MB), so **no data preparation is needed**.

1. Install **Python 3** if you don't have it: <https://www.python.org/downloads/>. On Windows, tick "Add python.exe to PATH" during install.
2. Download this repo (green **Code** button → **Download ZIP**, then unzip) or `git clone` it.
3. Start it:
   - **macOS:** double-click `start_mac.command`. The first time, macOS may block it: right-click → Open → Open.
   - **Windows:** double-click `start_windows.bat`.
   - **Any system, from a terminal:** `python3 -m pip install -r requirements.txt` then `python3 app.py`.

The first start installs the requirements (a minute or two). The app then opens at `http://localhost:8000` in your browser. To stop it, close the terminal window.

The model covers Kyiv and its immediate surroundings (about 50.20–50.62°N, 30.20–30.86°E, shown as a blue frame on the map). For another area, rebuild the model as described under *Rebuilding the model*.

In the app you can:

- **Place points** by clicking the map, or by typing or pasting coordinates (`50.4501, 30.5241`, `50°27'00"N 30°31'27"E`, `50°27.007'N 30°31.447'E`). Drag markers to move them. Names can be edited.
- **Set the antenna height for each point**, and tick **roof** if the antenna stands on a roof (height then counts from the roof top). Each point shows its ground elevation and the building height under it.
- **Set the radio**: frequency, handheld power (W), receiver sensitivity, antenna gain, losses.
- Press **Analyze**. Links are drawn in colour (green/yellow/orange/red; dashed = no line of sight). Click a link to get its **height profile** with buildings and the Fresnel zone.
- **Repeater search** (optional): ranked rooftop and ground sites, with each site's links to every point. Use `+` to add a candidate as a point.
- **Save/Load points** as CSV, and export the results as KML or CSV.

Points and settings are remembered in the browser. The app runs only on your computer (`127.0.0.1`); the browser just needs internet access for the map background.

## How it works

| Step | Script | What it does |
|---|---|---|
| 1 | `extract_buildings.py` | Reads building outlines and heights straight from an **OsmAnd `.obf` map** (OsmAnd's own binary format, decoded per its published spec `OBF.proto`). |
| 2 | `build_surface.py` | Merges **Copernicus GLO-30** terrain with the OSM buildings into a 10 m surface model. |
| 3 | `radio_los.py` | Runs the point-to-point analysis and the repeater search (command line). |
| 4 | `app.py` + `static/index.html` | The local map app built on the same analysis. |

### Building heights

The height for each building is chosen in this order:

1. The OSM `height` tag.
2. `building:levels` × 3 m + 1 m.
3. The median height of similar tagged buildings within 200 m (same class, footprint within ×2).
4. A default for the building type (house 6 m, apartment block 15 m, …).

This was tested on 5,000 Kyiv buildings whose height is known. Step 3 gives an average error of 6 m, against 8 m for type defaults alone. For towers above 40 m the error is still about 25 m, so a link blocked by a building with an *estimated* height is flagged in the output (`obstacle_height_from_osm_tag = estimated`).

About 22% of buildings in the Kyiv city area carry a height or floor count in OSM. Tall buildings are better covered than small ones.

Copernicus GLO-30 is officially a surface model, but over Kyiv it shows almost no building relief (1–6 m over 16–25-floor blocks). So it is used here as ground. A morphological opening removes the small remaining bumps.

### Radio model

- Earth curvature uses the standard 4/3 effective radius.
- **Fresnel clearance** is the fraction of the 1st Fresnel zone left free at the worst point. A value ≥ 0.6 counts as a clean path; a value < 0 means no line of sight.
- **Path loss** is free space plus knife-edge diffraction (ITU-R P.526, Deygout method, up to 7 edges).
- For **high site ↔ street-level** links (one end ≥ 30 m above ground, the other ≤ 10 m), the Okumura–Hata urban model is also computed and the lower loss is used. Pure diffraction ignores street reflections and is far too pessimistic in a city. The `model` column says which model was used.
- Link status: **good** ≥ 20 dB margin, **marginal** 10–20, **weak** 0–10, **no link** < 0.

Treat the results as planning guidance and confirm the best candidates with a field test.

## Rebuilding the model (optional)

Only needed for a different area or fresher OSM data. A model you build goes to `out/` and is used instead of the shipped `model/`.

### Data you need (not in the repo)

1. **OsmAnd map of Kyiv region**: `Ukraine_kyiv_europe.obf`, from OsmAnd → Download maps.
2. **Copernicus GLO-30 terrain tiles**: free, no account needed (AWS Open Data, maintained by Sinergise):
   - `https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N50_00_E030_00_DEM/Copernicus_DSM_COG_10_N50_00_E030_00_DEM.tif`
   - `https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N50_00_E029_00_DEM/Copernicus_DSM_COG_10_N50_00_E029_00_DEM.tif`

   Each tile is 1°×1° and named after its south-west corner. N50 E030 covers Kyiv city.

Put them in `data/`.

## Command-line use (without the app)

```bash
# 1. buildings from the OsmAnd map (~30 s)
python3 extract_buildings.py data/Ukraine_kyiv_europe.obf data/buildings.json

# 2. surface model (~30 s); area is BBOX in build_surface.py (Kyiv + margin)
python3 build_surface.py data/N50E029.tif data/N50E030.tif data/buildings.json out/kyiv

# 3. analysis
python3 radio_los.py points.csv --freq 446 --repeaters --out results/
```

### `points.csv`

```csv
name,lat,lon,height_m,on_roof
Base,50.4501,30.5241,12,yes
Team A,50.5080,30.4975,1.5,no
```

- `height_m` is the antenna height: above the ground (`on_roof = no`, the default) or above the roof at that spot (`on_roof = yes`).

### Main options

| Option | Default | Meaning |
|---|---|---|
| `--freq` | 446 | MHz |
| `--tx-dbm` | 37 | handheld power (37 dBm = 5 W) |
| `--sens-dbm` | -118 | receiver sensitivity |
| `--gain-tx`, `--gain-rx` | 0 | antenna gains, dBi |
| `--losses` | 3 | cable/body losses, dB |
| `--repeaters` | off | run the repeater search |
| `--rep-tx-dbm` | 44 | repeater power (44 dBm = 25 W) |
| `--rep-gain` | 3 | repeater antenna gain, dBi |
| `--rep-mast` | 3 | repeater antenna height above roof/ground, m |
| `--min-roof` | 25 | only roofs of buildings at least this tall (with a *tagged* height) are candidates |
| `--serve-margin` | 10 | margin (dB) that counts a point as served |
| `--search-pad` | 2000 | search this far (m) beyond the points' bounding box |

### Outputs (`results/`)

- `links.csv` is the full pair-by-pair table.
- `repeater_candidates.csv` lists the top 15 sites (at least 300 m apart), with the margin to each point.
- `radio_links.kml` / `radio_links.geojson` are for maps / QGIS.
- `radio_map.html` is an interactive map (needs internet for the OSM background).

The surface rasters in `model/` (or `out/`) (`kyiv_surface.tif`, `kyiv_ground.tif`, `kyiv_building_height.tif`, `kyiv_height_is_tagged.tif`) open directly in QGIS. `kyiv_surface.tif` also works as the elevation input for the QGIS **Visibility Analysis** plugin by Zoran Čučković, if you want viewsheds.

## Limits

- Trees are ignored, and so are indoor losses for radios inside buildings.
- Very tall entries can be towers or masts rather than roofs; candidates ≥ 190 m are labelled "check".
- OSM heights are only as good as the mappers made them. Estimated heights can be wrong by tens of metres for high-rises.
- To change the area, edit `BBOX` in `build_surface.py`. Above ~50 × 50 km at 10 m, memory use grows (~100 MB per raster).

## Data sources and licences

- **Buildings:** © [OpenStreetMap](https://www.openstreetmap.org/copyright) contributors, available under the Open Database Licence (ODbL). The rasters in `model/` are derived from OSM data and are shared under the ODbL.
- **Terrain:** produced using Copernicus WorldDEM™-30 © DLR e.V. 2010–2014 and © Airbus Defence and Space GmbH 2014–2018, provided under COPERNICUS by the European Union and ESA; all rights reserved. Accessed via the Registry of Open Data on AWS (managed by Sinergise).
- **OsmAnd file format:** decoded according to `OBF.proto` published in [osmandapp/OsmAnd-resources](https://github.com/osmandapp/OsmAnd-resources).
