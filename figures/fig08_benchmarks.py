"""Paper 1, Fig. 8 – nothing to calibrate against: (a) agreement of flood extent with the official FRMP EL11 hazard maps (critical success
index) for every run; (b) depth difference model - FRMP where both are wet (median and P10-P90); (c) Sentinel-1 new open water after the
11/01/2021 flood (passes of 14/01/2021) over SER-04, with the modelled January 2021 B extent and the FRMP T100 extent; (d) new open water in
the SER-04 model domain for every pass, with the delay after the flood peak.
Inputs: ../tables/raw/frmp_vs_model.csv, frmp_dmax_vs_model.csv, S1_SER04_comparison.csv; rasters in ../data: s1_newwater_2021-01-14_[AD].tif
(s1_vs_model.py), SER-04_p12_Jan2021B_G1_maxdepth.tif (HEC-RAS max depth, Jan 2021 B, C1), frmp_0100_S4.tif (FRMP EL11 T100 extent), model_perimeters.json.
usage: python fig08_benchmarks.py [data_dir]"""
import os, sys, csv, json
import numpy as np, rasterio
from rasterio.warp import reproject, Resampling
from shapely.geometry import Polygon
from pyproj import Transformer
from p1_style import plt, RAW, HERE, CASES, NAMES, CH_COL, save
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'data')
fv = list(csv.DictReader(open(os.path.join(RAW, 'frmp_vs_model.csv'))))
fd = list(csv.DictReader(open(os.path.join(RAW, 'frmp_dmax_vs_model.csv'))))
s1 = list(csv.DictReader(open(os.path.join(RAW, 'S1_SER04_comparison.csv'))))
fig = plt.figure(figsize=(13, 9.4)); gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1, 1.1])
MK = {'T50': 'o', 'T100': 's', 'T1000': 'D', 'T50CN3': '^'}
# (a)
ax = fig.add_subplot(gs[0, 0])
for i, c in enumerate(CASES):
    rr = [r for r in fv if r['case'] == c and r['run'] in MK]
    for j, r in enumerate(rr):
        ax.scatter(i + (j - len(rr) / 2) * .055, float(r['CSI']), marker=MK[r['run']], s=30, color=CH_COL[r['channel']], zorder=3)
ax.set_xticks(range(4)); ax.set_xticklabels([f'{c}\n{NAMES[c]}' for c in CASES]); ax.set_ylim(0, .65)
ax.set_ylabel('critical success index vs FRMP extent')
for t, m in MK.items(): ax.scatter([], [], marker=m, color='#555', label=t.replace('T50CN3', 'T50 CN III'))
for ch, col in CH_COL.items(): ax.scatter([], [], marker='o', color=col, label=ch)
ax.legend(ncol=4, fontsize=7.2, loc='upper left', frameon=False)
ax.set_title('(a) Extent agreement with the official FRMP maps (same T)', loc='left')
# (b)
ax = fig.add_subplot(gs[1, 0])
for i, c in enumerate(CASES):
    rr = [r for r in fd if r['case'] == c and r['run'] in ('T50', 'T100')]
    for j, r in enumerate(rr):
        x = i + (j - len(rr) / 2) * .09
        ax.plot([x, x], [float(r['dd_p10_m']), float(r['dd_p90_m'])], color=CH_COL[r['channel']], lw=1.3)
        ax.scatter(x, float(r['dd_median_m']), marker=MK[r['run']], color=CH_COL[r['channel']], s=28, zorder=3)
ax.axhline(0, color='#333', lw=.8); ax.axhspan(-.3, .3, color='#eeeeea', zorder=0)
ax.set_xticks(range(4)); ax.set_xticklabels([f'{c}\n{NAMES[c]}' for c in CASES]); ax.set_ylim(-1.7, 1.0)
ax.set_ylabel('max depth, model − FRMP (m), where both wet')
ax.set_title('(b) Depth agreement where both maps are wet (T50, T100)\nmarker = median, bar = P10–P90, band = ±0.3 m', loc='left')
# (c)
ax = fig.add_subplot(gs[0, 1])
with rasterio.open(os.path.join(DATA, 's1_newwater_2021-01-14_D.tif')) as r:
    W = r.read(1).astype(bool); tr = r.transform; crs = r.crs; H, Wd = W.shape; b = r.bounds
with rasterio.open(os.path.join(DATA, 's1_newwater_2021-01-14_A.tif')) as r: W |= r.read(1).astype(bool)
def grid(path, depth=True):
    with rasterio.open(path) as src:
        a = src.read(1); v = ((a >= .10) & (a != src.nodata)) if depth else (a > 0)
        o = np.zeros((H, Wd), 'uint8')
        reproject(v.astype('uint8'), o, src_transform=src.transform, src_crs=src.crs, dst_transform=tr, dst_crs=crs, resampling=Resampling.max)
    return o.astype(bool)
M = grid(os.path.join(DATA, 'SER-04_p12_Jan2021B_G1_maxdepth.tif'))
F = grid(os.path.join(DATA, 'frmp_0100_S4.tif'), depth=False)
ext = (b.left, b.right, b.bottom, b.top)
ax.imshow(np.ma.masked_where(~M, M), cmap='Blues', vmin=0, vmax=1.7, extent=ext, alpha=.6)
ax.contour(np.flipud(F.astype(float)), levels=[.5], colors='#e07b39', linewidths=.8, extent=ext, origin='lower')
from scipy import ndimage
yy, xx = np.nonzero(ndimage.binary_dilation(W, iterations=3)); X, Y = tr * (xx + .5, yy + .5)
ax.scatter(X, Y, s=1.2, c='crimson', marker='s', linewidths=0)
per = json.load(open(os.path.join(DATA, 'model_perimeters.json')))['SER-04']['pts']
t = Transformer.from_crs(2100, crs.to_epsg(), always_xy=True); px, py = zip(*[t.transform(x, y) for x, y in per])
ax.plot(px, py, 'k--', lw=.7); ax.set_xlim(min(px), max(px)); ax.set_ylim(min(py), max(py)); ax.set_aspect('equal')
ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
ax.plot([], [], 's', color='#7fb2dd', label='model Jan 2021 B (gauge rain), ≥ 0.10 m'); ax.plot([], [], '-', color='#e07b39', label='FRMP T100 extent')
ax.plot([], [], 's', color='crimson', label='Sentinel-1 new open water, 14/01/2021 (enlarged)'); ax.legend(loc='lower left', fontsize=7.2, frameon=True)
ax.set_title('(c) SER-04: radar 2.5–3 days after the 11/01/2021 flood', loc='left')
# (d)
ax = fig.add_subplot(gs[1, 1])
lab = {'2021-01-08': '8/1/21\n(+3 d, ev. A)', '2021-01-14': '14/1/21\n(+2.5–3 d, ev. B)', '2021-01-20': '20/1/21\n(+8.5 d)',
       '2026-01-24': '24/1/26\n(−2 d)', '2026-01-30': '30/1/26\n(+4 d)', '2026-02-05': '5/2/26\n(+10 d)'}
dates = list(lab); x = np.arange(len(dates))
for k, (tr_, col, nm) in enumerate([('D', '#2a78d6', 'descending (06:30 local)'), ('A', '#eb6834', 'ascending (18:24 local)')]):
    v = [float(next(r for r in s1 if r['date_track'] == f'{d}_{tr_}')['new_water_ha']) for d in dates]
    ax.bar(x + (k - .5) * .38, v, width=.36, color=col, label=nm)
    for xi, vi in zip(x + (k - .5) * .38, v): ax.text(xi, vi + .8, f'{vi:.1f}', ha='center', fontsize=7)
ax.set_xticks(x); ax.set_xticklabels([lab[d] for d in dates], fontsize=7.6); ax.set_ylabel('new open water in the model domain (ha)')
ax.annotate('scene-wide backscatter drop\n(~5 000 ha of the plain; not in the\nevening pass) – surface effect, not flood', xy=(2 - .19, 47), xytext=(3.1, 38),
            fontsize=7.4, arrowprops=dict(arrowstyle='->', color='#555'))
ax.text(.02, .62, 'model Jan 2021 B: 2 039 ha\nFRMP T100: 1 101 ha', transform=ax.transAxes, fontsize=7.6, color='#444')
ax.legend(loc='upper left', fontsize=7.5, frameon=False); ax.set_ylim(0, 55)
ax.set_title('(d) Sentinel-1 passes over SER-04 (delay after the flood peak)', loc='left')
fig.tight_layout(); save(fig, 'Fig08_benchmarks.png')
