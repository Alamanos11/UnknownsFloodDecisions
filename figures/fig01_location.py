"""Paper 1, Fig. 1 – study area: four catchments, model domains, streams, towns (a); the four model domains with the
reference flood extent (C1, T50, CN II, depth >= 0.10 m) and buildings (b-e).
usage: python fig01_location.py <02_CASES>   (needs the full project data: case DEMs, catchment gpkgs, result rasters, buildings – not in the repository)"""
import sys, os, json, glob
import numpy as np, rasterio, geopandas as gpd, pyogrio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from shapely.geometry import Polygon
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
ROOT = sys.argv[1]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Fig01_location.png')
CASES = [('SER-01', 'SER-01_AgioiAnargyroi', 'AgioiAnargyroi', 'Agioi Anargyroi', 'p01_T50_G2'),
         ('SER-02', 'SER-02_AgiaVarvara', 'AgiaVarvara', 'Agia Varvara', 'p01_T50_G1'),
         ('SER-03', 'SER-03_AgiosIoannis', 'AgiosIoannis', 'Agios Ioannis', 'p01_T50_G1'),
         ('SER-04', 'SER-04_Krousovitis', 'Krousovitis', 'Krousovitis', 'p01_T50_G1')]
COL = {'SER-01': '#6a3d9a', 'SER-02': '#b15928', 'SER-03': '#008080', 'SER-04': '#d6457f'}   # = p1_style.CASE_COL
per = json.load(open(os.path.join(ROOT, 'SER-P1_paper_decisions', 'data', 'model_perimeters.json')))
# mosaic hillshade 25 m
x0, x1, y0, y1, px = 444000, 474000, 4541000, 4584000, 25
W, H = int((x1 - x0) / px), int((y1 - y0) / px); tr = from_origin(x0, y1, px, px)
Z = np.full((H, W), np.nan, 'float32')
for c, d, n, *_ in CASES:
    with rasterio.open(os.path.join(ROOT, d, '00_data', 'dem', f'{n}_dem5m_filled.tif')) as r:
        a = np.full((H, W), np.nan, 'float32')
        reproject(r.read(1).astype('float32'), a, src_transform=r.transform, src_crs=r.crs, dst_transform=tr, dst_crs=r.crs,
                  src_nodata=r.nodata, dst_nodata=np.nan, resampling=Resampling.average)
        Z = np.where(np.isnan(Z), a, Z)
ls = LightSource(azdeg=315, altdeg=45)
hs = ls.hillshade(np.nan_to_num(Z, nan=np.nanmin(Z)), vert_exag=2, dx=px, dy=px)
fig = plt.figure(figsize=(13, 9.2))
ax = fig.add_axes([0.02, 0.05, 0.42, 0.9])
ext = (x0, x1, y0, y1)
ax.imshow(np.ma.masked_where(np.isnan(Z), hs), cmap='gray', extent=ext, vmin=0, vmax=1.2, alpha=.9)
ax.imshow(np.ma.masked_where(np.isnan(Z), Z), cmap='terrain', extent=ext, alpha=.25, vmin=0, vmax=1800)
for c, d, n, name, _ in CASES:
    g = os.path.join(ROOT, d, '01_gis', 'outputs', f'{n}_catchment_5m.gpkg')
    cat = gpd.read_file(g, layer='catchment'); st = gpd.read_file(g, layer='streams_0p1km2'); ou = gpd.read_file(g, layer='outlet')
    cat.boundary.plot(ax=ax, color=COL[c], lw=1.6)
    st.plot(ax=ax, color='#3b6fb6', lw=.35, alpha=.7)
    ou.plot(ax=ax, color='k', markersize=14, zorder=5)
    gpd.GeoSeries([Polygon(per[c]['pts'])]).boundary.plot(ax=ax, color='k', lw=.9, ls='--')
    cx, cy = cat.geometry.iloc[0].representative_point().coords[0]
    ax.text(cx, cy, f'{c}\n{name}', ha='center', va='center', fontsize=8.5, color='k',
            bbox=dict(boxstyle='round,pad=.25', fc='white', ec=COL[c], lw=1, alpha=.85))
for t, xy, off in [('Serres', (461300, 4549300), (-4200, -1600)), ('Sidirokastro', (447600, 4566300), (-3300, 900))]:
    ax.plot(*xy, 's', color='#333', ms=5); ax.text(xy[0] + off[0], xy[1] + off[1], t, fontsize=8.5, style='italic', bbox=dict(fc='white', ec='none', alpha=.7, pad=1))
ax.plot([], [], color='k', ls='--', lw=.9, label='2D model domain'); ax.plot([], [], color='#555', lw=1.6, label='catchment')
ax.plot([], [], 'o', color='k', ms=4, label='outlet (hydrology)')
ax.legend(loc='lower left', fontsize=8, frameon=True)
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(''); ax.set_ylabel('')
ax.plot([x0 + 1500, x0 + 6500], [y0 + 4200] * 2, 'k-', lw=3); ax.text(x0 + 4000, y0 + 4700, '5 km', ha='center', fontsize=8)
ax.annotate('N', xy=(x1 - 1500, y1 - 1500), xytext=(x1 - 1500, y1 - 4000), ha='center', fontsize=10,
            arrowprops=dict(arrowstyle='-|>', color='k'))
ax.set_title('(a) Four torrents, Serres regional unit (EL11), EGSA87', fontsize=10, loc='left')
# panels b-e
bld = gpd.read_file(os.path.join(ROOT, '_data_buildings', 'buildings_merged_serres.gpkg'))
pos = [[0.47, 0.52, 0.25, 0.43], [0.74, 0.52, 0.25, 0.43], [0.47, 0.05, 0.25, 0.43], [0.74, 0.05, 0.25, 0.43]]
for (c, d, n, name, plan), p, lab in zip(CASES, pos, 'bcde'):
    a2 = fig.add_axes(p); poly = Polygon(per[c]['pts']); bx = poly.bounds
    f = glob.glob(os.path.join(ROOT, d, '03_hecras', 'results', plan, '*maxdepth.tif'))[0]
    with rasterio.open(f) as r:
        D = r.read(1); D = np.where((D >= 0.10) & (D != r.nodata), 1, np.nan); e = r.bounds
    a2.imshow(D, extent=(e.left, e.right, e.bottom, e.top), cmap='Blues', vmin=0, vmax=1.3, interpolation='nearest')
    b = bld.cx[bx[0]:bx[2], bx[1]:bx[3]]
    b.plot(ax=a2, color='#555', lw=0, markersize=1)
    gpd.GeoSeries([poly]).boundary.plot(ax=a2, color='k', lw=.8, ls='--')
    a2.set_xlim(bx[0], bx[2]); a2.set_ylim(bx[1], bx[3]); a2.set_aspect('equal'); a2.set_xticks([]); a2.set_yticks([]); a2.set_xlabel(''); a2.set_ylabel('')
    for s in a2.spines.values(): s.set_edgecolor(COL[c]); s.set_linewidth(1.5)
    L = 1000 if c != 'SER-04' else 2000
    a2.plot([bx[0] + 0.05 * (bx[2] - bx[0]), bx[0] + 0.05 * (bx[2] - bx[0]) + L], [bx[1] + 0.03 * (bx[3] - bx[1])] * 2, 'k-', lw=2.5)
    a2.text(bx[0] + 0.05 * (bx[2] - bx[0]) + L / 2, bx[1] + 0.045 * (bx[3] - bx[1]), f'{L // 1000} km', ha='center', va='bottom', fontsize=7.5, bbox=dict(fc='white', ec='none', alpha=.7, pad=1))
    a2.set_title(f'({lab}) {c} {name} – domain {poly.area / 1e6:.1f} km²', fontsize=9, loc='left')
a2.plot([], [], 's', color='#6fa8dc', label='flood extent C1 T50 (≥ 0.10 m)'); a2.plot([], [], 's', color='#555', ms=4, label='buildings')
a2.legend(loc='upper right', fontsize=7.5, frameon=True)
fig.savefig(OUT, dpi=300)
print('ok')
