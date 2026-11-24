# Which unknowns matter for which flood decision?

Supporting material for the paper *"Which unknowns matter for which flood decision: Decision-relevant uncertainty in survey-free 2D
modelling of ungauged torrents"*

by Angelos Alamanos

Four ungauged torrents in Serres (Greece) were modelled with HEC-HMS and HEC-RAS 2D over a factorial of unknowns (channel geometry,
antecedent soil moisture, design return period) and management levers (maintenance states), and each run was scored against four
decisions: D1 which stream/reach to study first, D2 whether maintenance is worth it, D3 at what discharge damage starts, D4 where the
flood line falls. A two-way variance decomposition shows which unknown matters for which decision.

The HEC-HMS / HEC-RAS projects are not included (size). This repository holds the figures, the tables and the result files behind them,
and the scripts that produced the analyses.

## Contents
| Folder | What |
|---|---|
| `figures/` | Figures 1–8 (PNG, 300 dpi) and the script that draws each (`fig0N_*.py`, shared style `p1_style.py`) |
| `tables/` | Tables 1–5 (CSV; `make_tables.py`, `triage_costs.py`) and `raw/` – the result tables of every analysis (CSV) |
| `data/` | small rasters used by Figs 6 and 8: number of plausible states flooding each pixel (`<case>_T50_state_count.tif`), Sentinel-1 new-water masks of 14/01/2021, the SER-04 January 2021 B max-depth raster, the FRMP T100 extent (SER-04 area), model perimeters (EPSG:2100) |
| `scripts/core/` | analysis: `p1_common.py` (case configuration), `d1_priorities.py` (D1), `scen_post.py` (D2), `d3.py`, `d3_summary.py` (D3), `d4_floodlines.py` (D4), `anova_decisions.py` (variance decomposition) |
| `scripts/frmp/` | reading the official FRMP EL11 maps from the national geoportal WMS (`wms_dmax_decode.js`, run in a browser) and the comparison with the model (`frmp_compare.py`, `dmax_compare.py`) |
| `scripts/sentinel1/` | Sentinel-1 change detection over SER-04 (`s1_browser_fetch.js`, run in the Copernicus Browser; `s1_masks.json`; `s1_vs_model.py`) |

## Reproducing
- Python 3.10+ with numpy, scipy, matplotlib, rasterio, geopandas, shapely, pyproj, pandas.
- **From this repository alone:** `cd figures; python fig02_design.py` … `python fig08_benchmarks.py` (Figs 2–8) and
  `cd tables; python make_tables.py` redraw the figures and tables from `tables/raw/` and `data/`.
  Fig. 1 needs the full project data (DEMs, catchments, buildings) and is provided as an image.
- **The analyses** (`scripts/`) need the HEC-RAS result rasters and building footprints of the project: set the environment variable
  `P1_ROOT` to the project's `02_CASES` folder (and, for D2, `SCEN_TERRAIN` and `BLD_DIR`). They document the method exactly as run.

## Software and data
HEC-HMS 4.12 and HEC-RAS 6.7 (US Army Corps of Engineers, HEC). Design rainfall: national IDF parameters (NTUA, 2023).
Buildings: OpenStreetMap (© OpenStreetMap contributors, ODbL) and Microsoft Global ML Building Footprints (ODbL).
Benchmarks: FRMP EL11 hazard maps (Ministry of Environment and Energy, floods.ypeka.gr; derived extents and depths);
Sentinel-1 IW GRD via the Copernicus Data Space Ecosystem (contains modified Copernicus Sentinel data 2020–2026).
