"""Paper 1, Fig. 6 – D4 where the line falls: (a-d) number of the 14 plausible T50 states (channel C0/C1/C+, AMC II/III, maintenance)
that flood each pixel (>= 0.10 m) – 14 = always, 1-13 = the band where the legal line may fall; (e) shift of the T50 line from the
reference (C1, CN II) for one unknown at a time: median and P90 of the absolute displacement.
Inputs: ../data/<case>_T50_state_count.tif (from d4_floodlines.py), ../data/model_perimeters.json;
../tables/raw/D4_line_displacement.csv, D4_band_summary.csv.   usage: python fig06_d4_band.py [data_dir]"""
import os, sys, csv, json
import numpy as np, rasterio
from shapely.geometry import Polygon
from matplotlib.colors import ListedColormap, BoundaryNorm
from p1_style import plt, RAW, HERE, CASES, NAMES, CASE_COL, save
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'data')
per = json.load(open(os.path.join(DATA, 'model_perimeters.json')))
band = {r['case']: r for r in csv.DictReader(open(os.path.join(RAW, 'D4_band_summary.csv')))}
disp = list(csv.DictReader(open(os.path.join(RAW, 'D4_line_displacement.csv'))))
cols = ['#fde0c5', '#f9b98a', '#f08a5d', '#d45a3a', '#a33224', '#2a5aa0']   # 1-2, 3-5, 6-8, 9-11, 12-13 states, 14 = always
cmap = ListedColormap(cols); norm = BoundaryNorm([.5, 2.5, 5.5, 8.5, 11.5, 13.5, 14.5], cmap.N)
fig = plt.figure(figsize=(13, 9.6))
pos = [[.03, .40, .22, .55], [.27, .40, .22, .55], [.51, .40, .22, .55], [.75, .40, .22, .55]]
for (c, p), lab in zip(zip(CASES, pos), 'abcd'):
    ax = fig.add_axes(p)
    with rasterio.open(os.path.join(DATA, f'{c}_T50_state_count.tif')) as r:
        A = r.read(1).astype(float); e = r.bounds; nd = r.nodata
    A[(A <= 0) | (A > 14) | (A == (nd if nd is not None else -1))] = np.nan
    ax.imshow(A, extent=(e.left, e.right, e.bottom, e.top), cmap=cmap, norm=norm, interpolation='nearest')
    poly = Polygon(per[c]['pts']); x, y = poly.exterior.xy; ax.plot(x, y, color='#333', lw=.7, ls='--')
    bx = poly.bounds; ax.set_xlim(bx[0], bx[2]); ax.set_ylim(bx[1], bx[3]); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for s in ax.spines.values(): s.set_visible(True); s.set_edgecolor(CASE_COL[c])
    b = band[c]
    ax.set_title(f'({lab}) {c} {NAMES[c]}\nband {float(b["A_sometimes_ha"]):.0f} ha, {b["bld_sometimes"]} buildings', loc='left', fontsize=9)
    L = 1000 if c != 'SER-04' else 2000; x0 = bx[0] + .05 * (bx[2] - bx[0]); y0 = bx[1] + .03 * (bx[3] - bx[1])
    ax.plot([x0, x0 + L], [y0, y0], 'k-', lw=2.5); ax.text(x0 + L / 2, y0 + .015 * (bx[3] - bx[1]), f'{L // 1000} km', ha='center', fontsize=7.5)
cax = fig.add_axes([.30, .355, .40, .015])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal', ticks=[1.5, 4, 7, 10, 12.5, 14])
cb.ax.set_xticklabels(['1–2', '3–5', '6–8', '9–11', '12–13', '14 (always)'], fontsize=8)
cb.set_label('number of the 14 plausible T50 states that flood the pixel (depth ≥ 0.10 m)', fontsize=8.5)
ax = fig.add_axes([.07, .05, .88, .24])
UNK = [('C1 T100', 'design standard T100'), ('C1 T50 CN III', 'soil moisture AMC III'), ('C0 T50', 'channel C0 (DEM only)'), ('C+ T50', 'channel C+ (larger)')]
UC = ['#9aa0a6', '#eb6834', '#e07b39', '#1baf7a']
w = .2
for k, ((sc, lab), col) in enumerate(zip(UNK, ['#9aa0a6', '#7b61c4', '#e07b39', '#1baf7a'])):
    med, p90 = [], []
    for c in CASES:
        r = next(x for x in disp if x['case'] == c and x['scenario'] == sc)
        med.append(float(r['median_abs_disp_m'])); p90.append(float(r['p90_abs_disp_m']))
    xx = np.arange(4) + (k - 1.5) * w
    ax.bar(xx, p90, width=w * .9, color=col, alpha=.35, edgecolor=col)
    ax.bar(xx, med, width=w * .9, color=col, label=lab)
    for xi, m, p9 in zip(xx, med, p90): ax.text(xi, p9 + 8, f'{p9:.0f}', ha='center', fontsize=6.8, color='#444')
ax.set_xticks(range(4)); ax.set_xticklabels([f'{c} {NAMES[c]}' for c in CASES])
ax.set_ylabel('line shift from reference (m)'); ax.set_ylim(0, 980)
ax.legend(ncol=2, loc='upper left', frameon=False, fontsize=8)
ax.text(.42, .97, 'solid = median, light = P90 (value printed)', transform=ax.transAxes, ha='left', va='top', fontsize=7.5, color='#555')
ax.set_title('(e) How far the T50 line moves when one unknown changes (reference: C1, CN II, T50)', loc='left', fontsize=9.5)
save(fig, 'Fig06_D4_band.png')
