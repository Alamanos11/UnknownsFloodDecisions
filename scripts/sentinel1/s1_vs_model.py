"""Sentinel-1 check of SER-04 (Krousovitis) flood extents.

Input  : s1_masks.json  (bit-packed masks on an EPSG:3857 grid, made in the Copernicus Browser
         session with the page's own Sentinel Hub layer; see METHODS.md section 4)
         model maxdepth tifs (EPSG:2100, 5 m), FRMP EL11 extent tifs, model_perimeters.json
Output : S1_SER04_comparison.csv, S1_SER04_clusters.csv, fig_S1_SER04.png, s1_newwater_<date>.tif

usage: python s1_vs_model.py <CASES_ROOT> <OUT_DIR>
"""
import sys, json, gzip, base64, os
import numpy as np, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds
from rasterio.features import rasterize
from scipy import ndimage
from shapely.geometry import Polygon, mapping
from shapely.ops import transform as shp_tr
from pyproj import Transformer
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 's1_masks.json')))
m = d['meta']; W, H = m['W'], m['H']
TR = from_bounds(m['x0'], m['y0'], m['x1'], m['y1'], W, H)
PXHA = None  # ground area per pixel computed below

def unpack(s):
    b = gzip.decompress(base64.b64decode(s))
    return np.unpackbits(np.frombuffer(b, np.uint8))[:W * H].reshape(H, W).astype(bool)

masks = {k: unpack(v) for k, v in d['masks'].items()}

# ground pixel area (Mercator scale at each row)
R = 6378137.0
ys = m['y1'] - (np.arange(H) + 0.5) * (m['y1'] - m['y0']) / H
lat = np.degrees(2 * np.arctan(np.exp(ys / R)) - np.pi / 2)
pxm = (m['x1'] - m['x0']) / W
area_ha = ((pxm * np.cos(np.radians(lat))) ** 2 / 1e4)[:, None] * np.ones((1, W))

# model perimeter -> S1 grid
per = json.load(open(os.path.join(ROOT, 'SER-P1_paper_decisions', 'data', 'model_perimeters.json')))['SER-04']['pts']
t = Transformer.from_crs(2100, 3857, always_xy=True)
poly = shp_tr(lambda x, y: t.transform(x, y), Polygon(per))
inside = rasterize([(mapping(poly), 1)], out_shape=(H, W), transform=TR, fill=0).astype(bool)

def to_grid(path, wet=0.10, is_depth=True):
    with rasterio.open(path) as src:
        a = src.read(1).astype('float32')
        nd = src.nodata
        if is_depth:
            v = np.where((a >= wet) & (a != nd), 1, 0).astype('uint8')
        else:
            v = (a > 0).astype('uint8')
        out = np.zeros((H, W), 'uint8')
        reproject(v, out, src_transform=src.transform, src_crs=src.crs, dst_transform=TR,
                  dst_crs='EPSG:3857', resampling=Resampling.max)
    return out.astype(bool)

res = os.path.join(ROOT, 'SER-04_Krousovitis', '03_hecras', 'results')
frmp = os.path.join(ROOT, '_data_frmp', 'frmp_el11_extents')
refs = {
    'C1_T50': to_grid(os.path.join(res, 'p01_T50_G1', 'maxdepth.tif')),
    'C1_T100': to_grid(os.path.join(res, 'p02_T100_G1', 'maxdepth.tif')),
    'C1_T1000': to_grid(os.path.join(res, 'p03_T1000_G1', 'maxdepth.tif')),
    'C0_T50': to_grid(os.path.join(res, 'p06_T50_G0', 'maxdepth.tif')),
    'Cp_T50': to_grid(os.path.join(res, 'p07_T50_G2', 'maxdepth.tif')),
    'Jan2021A_gauge': to_grid(os.path.join(res, 'p05_Jan2021A_G1', 'maxdepth.tif')),
    'Jan2021A_cerra': to_grid(os.path.join(res, 'p13_Jan2021Ac_G1', 'maxdepth.tif')),
    'Jan2021B_gauge': to_grid(os.path.join(res, 'p12_Jan2021B_G1', 'maxdepth.tif')),
    'FRMP_T50': to_grid(os.path.join(frmp, 'frmp_0050_S4.tif'), is_depth=False),
    'FRMP_T100': to_grid(os.path.join(frmp, 'frmp_0100_S4.tif'), is_depth=False),
    'FRMP_T1000': to_grid(os.path.join(frmp, 'frmp_1000_S4.tif'), is_depth=False),
}
# 1-pixel tolerance (13 m) for edge registration
tol = {k: ndimage.binary_dilation(v, iterations=1) for k, v in refs.items()}

def clean(mk):
    # drop isolated speckle: keep connected groups >= 5 px (~0.05 ha)
    lab, n = ndimage.label(mk, structure=np.ones((3, 3)))
    if n == 0:
        return mk
    sz = ndimage.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1))
    keep = np.zeros(n + 1, bool); keep[1:] = sz >= 5
    return keep[lab]

rows = []
event_dates = [k for k in masks if not k.startswith('refdark')]
for k in sorted(event_dates):
    mk = clean(masks[k]) & inside & ~ndimage.binary_dilation(masks['refdark_D21'], iterations=1)
    a = float(area_ha[mk].sum())
    r = {'date_track': k, 'new_water_ha': round(a, 2)}
    for rn, rv in tol.items():
        r['in_' + rn + '_pct'] = round(100 * area_ha[mk & rv].sum() / a, 1) if a > 0 else np.nan
    r['modelJan2021B_detected_pct'] = round(100 * area_ha[mk & refs['Jan2021B_gauge'] & inside].sum() /
                                            area_ha[refs['Jan2021B_gauge'] & inside].sum(), 2)
    rows.append(r)
    with rasterio.open(os.path.join(OUT, f's1_newwater_{k}.tif'), 'w', driver='GTiff', height=H, width=W,
                       count=1, dtype='uint8', crs='EPSG:3857', transform=TR, compress='deflate') as dst:
        dst.write(mk.astype('uint8'), 1)

import csv
with open(os.path.join(OUT, 'S1_SER04_comparison.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# model areas inside perimeter for context
ctx = {k: round(float(area_ha[v & inside].sum()), 1) for k, v in refs.items()}
json.dump({'model_and_frmp_area_ha_on_S1_grid': ctx, 'refdark_ha': round(float(area_ha[masks['refdark_D21'] & inside].sum()), 1)},
          open(os.path.join(OUT, 'S1_SER04_context.json'), 'w'), indent=1)

# clusters of the largest signal (20 Jan descending) and of 14 Jan
inv = Transformer.from_crs(3857, 2100, always_xy=True)
crow = []
for k in ['2021-01-14_D', '2021-01-14_A', '2021-01-20_D', '2021-01-20_A', '2026-01-30_D']:
    mk = clean(masks[k]) & inside & ~ndimage.binary_dilation(masks['refdark_D21'], iterations=1)
    lab, n = ndimage.label(mk, structure=np.ones((3, 3)))
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        sel = lab[sl] == i
        a = float(area_ha[sl][sel].sum())
        if a < 0.5:
            continue
        cy, cx = ndimage.center_of_mass(lab == i)
        X, Y = TR * (cx + 0.5, cy + 0.5)
        E, N = inv.transform(X, Y)
        crow.append({'date_track': k, 'cluster': i, 'area_ha': round(a, 2), 'E_2100': round(E), 'N_2100': round(N),
                     'in_Jan2021B_model': bool((tol['Jan2021B_gauge'][sl][sel]).mean() > 0.5),
                     'in_C1_T100': bool((tol['C1_T100'][sl][sel]).mean() > 0.5),
                     'in_FRMP_T100': bool((tol['FRMP_T100'][sl][sel]).mean() > 0.5)})
with open(os.path.join(OUT, 'S1_SER04_clusters.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(crow[0].keys())); w.writeheader(); w.writerows(sorted(crow, key=lambda r: (r['date_track'], -r['area_ha'])))

# figure
fig, axs = plt.subplots(1, 3, figsize=(16, 6.4), sharex=True, sharey=True)
ext = (m['x0'], m['x1'], m['y0'], m['y1'])
for ax, (k, ttl) in zip(axs, [('2021-01-14_D', '14/1/2021 06:31 local (desc.)\n~2.5 days after the 11/1 peak'),
                               ('2021-01-20_D', '20/1/2021 06:30 local (desc.)\n~8.5 days after the peak'),
                               ('2026-01-30_D', '30/1/2026 06:30 local (desc.)\n4 days after the 26/1 flood')]):
    ax.imshow(np.where(inside, 0.93, 1.0), cmap='gray', vmin=0, vmax=1, extent=ext)
    ax.imshow(np.ma.masked_where(~refs['Jan2021B_gauge'], refs['Jan2021B_gauge']), cmap='Blues', vmin=0, vmax=1.6, alpha=.55, extent=ext)
    ax.contour(np.flipud(refs['FRMP_T100'].astype(float)), levels=[.5], colors='darkorange', linewidths=.6,
               extent=ext, origin='lower')
    mk = clean(masks[k]) & inside & ~ndimage.binary_dilation(masks['refdark_D21'], iterations=1)
    yy, xx = np.nonzero(ndimage.binary_dilation(mk, iterations=2))
    X, Y = TR * (xx + 0.5, yy + 0.5)
    ax.scatter(X, Y, s=.4, c='crimson', marker='s', linewidths=0)
    ax.set_title(ttl, fontsize=10); ax.set_xticks([]); ax.set_yticks([])
axs[0].plot([], [], 's', color='#6fa8dc', label='model Jan 2021 B (gauge, C1), depth ≥ 0.1 m')
axs[0].plot([], [], '-', color='darkorange', label='FRMP EL11 T100 extent')
axs[0].plot([], [], 's', color='crimson', label='Sentinel-1 new open water (VV ≤ −18 dB, drop ≥ 3 dB)')
axs[0].legend(loc='lower left', fontsize=8, frameon=True)
fig.suptitle('SER-04 Krousovitis – Sentinel-1 VV change detection vs model and FRMP (grey = model domain)', fontsize=11, y=0.995)
fig.tight_layout(rect=(0,0,1,0.95))
fig.savefig(os.path.join(OUT, 'fig_S1_SER04.png'), dpi=170)
print(json.dumps(ctx)); [print(r) for r in rows]
