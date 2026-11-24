"""D4 - where does the flood line fall, and how far does it move with each unknown?

Flood extent = max depth >= 0.10 m inside the model domain, cleaned as the delivered flood lines were
(binary opening 1 cell, holes filled, patches < 400 m2 dropped, closing with a 4 m radius).
Flood line = extent pixels with a dry 4-neighbour; pixels within 3 cells of the domain edge are left out
(there the extent is cut by the model boundary, not by the terrain). Line length = length of the vectorised
extent boundary (simplified 1 m) inside the domain shrunk by 3 cells.
Buildings 'inside the line' = footprint touches the cleaned extent.
Reference line: C1, T50, CN II. For every other scenario (and the FRMP EL11 T50/T100 extents):
  displacement at each reference-line pixel = distance to the nearest line pixel of the scenario,
    signed + where the scenario floods beyond the reference line (pixel lies inside the scenario extent), - otherwise;
  reported as median and 90th percentile of |displacement| and the share of the reference line moving > 10 / 25 / 50 m;
  band = area flooded in one map but not the other; mean band width = band area / reference line length
    (line length ~ number of line pixels x cell size);
  buildings changing side = buildings reached (>= 0.10 m on the footprint) in one map but not in the other.
Ensemble of plausible T50 states (C0 / C1 / C+ x CN II / III, and C1 / C+ x 4 maintenance states = 14 maps):
  number of states that flood each pixel (0-14) -> always / sometimes / never; uncertain band per reach
  (area flooded in 1-13 states / reach length) and buildings in the band.
Per reach (reach of the nearest centreline station): median |displacement| for each single-factor change.
Writes SER-P1_paper_decisions/D4_floodlines/: D4_line_displacement.csv, D4_reach_displacement.csv, D4_band_summary.csv,
D4_band_reaches.csv, <case>_T50_state_count.tif, fig_D4_band_<case>.png, fig_D4_displacement.png.
usage: P1_ROOT=<02_CASES> python3 d4_floodlines.py"""
import os, csv, numpy as np
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.ops import unary_union
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from p1_common import CASES, REF, MAINT, Case, skey, P1, WET

OUT = os.path.join(P1, 'D4_floodlines'); os.makedirs(OUT, exist_ok=True)
ENS = [(ch, 'T50', cn, 'base') for ch in ('C0', 'C1', 'C+') for cn in ('II', 'III')] + [(ch, 'T50', 'II', m) for ch in ('C1', 'C+') for m in MAINT]
FACT = [('C0 (no channel)', ('C0', 'T50', 'II', 'base')), ('C+ (larger channel)', ('C+', 'T50', 'II', 'base')),
        ('CN III (wet soil)', ('C1', 'T50', 'III', 'base')), ('T100 (design standard)', ('C1', 'T100', 'II', 'base'))] + \
       [(f'{m} (maintenance)', ('C1', 'T50', 'II', m)) for m in MAINT]
L1, L2, L3, L4 = [], [], [], []
def clean(E, k):
    E = ndi.binary_opening(E, iterations=1); E = ndi.binary_fill_holes(E)
    lab, n = ndi.label(E)
    if n:
        sz = ndi.sum(E, lab, range(1, n + 1)) * k.px * k.px; E = np.isin(lab, 1 + np.where(sz >= 400)[0])
    r = max(1, int(round(4 / k.px))); st = ndi.iterate_structure(ndi.generate_binary_structure(2, 1), r)
    return ndi.binary_closing(np.pad(E, r), structure=st)[r:-r, r:-r] & k.dom
def line_km(E, k, inner_poly):
    g = unary_union([shape(p) for p, v in shapes(E.astype('uint8'), mask=E, transform=k.tr) if v == 1])
    return g.simplify(1.0).boundary.intersection(inner_poly).length / 1000
def line_of(E, inner):
    return E & ~ndi.binary_erosion(E, structure=ndi.generate_binary_structure(2, 1)) & inner
for c in CASES:
    k = Case(c); px = k.px; b, _ = k.buildings(); nr = len(k.edges) - 1; rt = k.reach_table()
    inner = ndi.binary_erosion(k.dom, iterations=3)
    inner_poly = k.poly.buffer(-3 * px)
    Dref = k.depth(REF); Eref = clean(Dref >= WET, k); Lref = line_of(Eref, inner); Bref = k.bld_depth(Eref.astype('float32')) > 0
    lref_km = line_km(Eref, k, inner_poly); rref = k.reach[Lref]
    maps = {skey(key): clean(k.depth(key) >= WET, k) for key in k.cf['plans'] if key != REF}
    maps['FRMP EL11 T50'] = clean(k.frmp('0050'), k); maps['FRMP EL11 T100'] = clean(k.frmp('0100'), k)
    disp = {}
    for name, E in maps.items():
        Ls = line_of(E, inner)
        dt = ndi.distance_transform_edt(~Ls) * px
        d = dt[Lref]; sgn = np.where(E[Lref] & ~Ls[Lref], 1, -1); disp[name] = d * sgn; a = np.abs(d)
        Bs = k.bld_depth(E.astype('float32')) > 0
        L1.append(dict(case=c, scenario=name, ref_line_km=round(lref_km, 2), A_ref_ha=round(Eref.sum() * k.pxa, 1), A_ha=round(E.sum() * k.pxa, 1),
                       median_abs_disp_m=round(float(np.median(a)), 1), p90_abs_disp_m=round(float(np.percentile(a, 90)), 1),
                       share_gt10m=round(float((a > 10).mean()), 3), share_gt25m=round(float((a > 25).mean()), 3), share_gt50m=round(float((a > 50).mean()), 3),
                       share_outward=round(float(((d > 0) & (sgn > 0)).mean()), 3),
                       band_ha=round(float((E ^ Eref).sum() * k.pxa), 1), mean_band_width_m=round(float((E ^ Eref).sum() * px * px / (lref_km * 1000)), 1),
                       bld_newly_in=int((Bs & ~Bref).sum()), bld_newly_out=int((~Bs & Bref).sum()), bld_ref=int(Bref.sum())))
    for j in range(nr):
        m = rref == j
        if m.sum() < 5: continue
        row = dict(case=c, reach=rt[j]['reach'], ch_from=rt[j]['ch_from'], ch_to=rt[j]['ch_to'], ref_line_pixels=int(m.sum()))
        for lab, key in FACT:
            if skey(key) in disp: row[f'med_abs_m {lab}'] = round(float(np.median(np.abs(disp[skey(key)][m]))), 1)
        row['med_abs_m FRMP T50'] = round(float(np.median(np.abs(disp['FRMP EL11 T50'][m]))), 1)
        L2.append(row)
    # ensemble of plausible T50 states
    ens = [key for key in ENS if key in k.cf['plans']]; N = np.zeros(k.sh, 'uint8'); BN = np.zeros(len(b), 'uint8')
    for key in ens:
        E = clean(k.depth(key) >= WET, k); N += E; BN += (k.bld_depth(E.astype('float32')) > 0)
    n = len(ens); always = (N == n) & k.dom; some = (N > 0) & (N < n) & k.dom
    k.write(os.path.join(OUT, f'{c}_T50_state_count.tif'), np.where(k.dom, N, 255), 'uint8', 255)
    L3.append(dict(case=c, states=n, A_ref_ha=round(Eref.sum() * k.pxa, 1), A_always_ha=round(always.sum() * k.pxa, 1),
                   A_sometimes_ha=round(some.sum() * k.pxa, 1), A_envelope_ha=round((N > 0).sum() * k.pxa, 1),
                   bld_always=int((BN == n).sum()), bld_sometimes=int(((BN > 0) & (BN < n)).sum()), bld_ref=int(Bref.sum()),
                   bld_in_domain=len(b)))
    for j in range(nr):
        L = rt[j]['length_m'];
        if L <= 0: continue
        sel = k.reach == j; bs = b['reach_i'].values == j
        L4.append(dict(case=c, reach=rt[j]['reach'], ch_from=rt[j]['ch_from'], ch_to=rt[j]['ch_to'], length_m=L,
                       always_ha=round((always & sel).sum() * k.pxa, 2), sometimes_ha=round((some & sel).sum() * k.pxa, 2),
                       band_width_m=round((some & sel).sum() * px * px / L, 1),
                       bld_always=int((BN[bs] == n).sum()), bld_sometimes=int(((BN[bs] > 0) & (BN[bs] < n)).sum())))
    # map of the state count with the reference line
    fig, ax = plt.subplots(figsize=(8, 8 * k.sh[0] / k.sh[1] if k.sh[1] > k.sh[0] else 8))
    ext = [k.tr.c, k.tr.c + k.tr.a * k.sh[1], k.tr.f + k.tr.e * k.sh[0], k.tr.f]
    cm = ListedColormap(['#f2f2f2'] + list(plt.cm.YlOrRd(np.linspace(0.15, 0.85, n - 1))) + ['#08306b'])
    A = np.where(k.dom, N, np.nan).astype(float)
    im = ax.imshow(A, extent=ext, cmap=cm, norm=BoundaryNorm(np.arange(-0.5, n + 1.5), n + 1), interpolation='nearest')
    ax.contour(np.flipud(Eref.astype(float)), levels=[0.5], extent=[ext[0], ext[1], ext[2], ext[3]], colors='k', linewidths=0.5, origin='lower')
    bx, by = b.representative_point().x, b.representative_point().y
    ax.scatter(bx[(BN > 0) & (BN < n)], by[(BN > 0) & (BN < n)], s=1.5, c='#6a3d9a', label='buildings in the uncertain band')
    ax.plot(*k.poly.exterior.xy, 'k-', lw=0.5); ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    cb = plt.colorbar(im, ax=ax, fraction=0.035, ticks=range(0, n + 1, 2)); cb.set_label(f'number of plausible T50 states flooding the cell (of {n})')
    ax.set_title(f"{c} {k.cf['name']}: T50 flood extent over {n} plausible states; black = reference line (C1, CN II)", fontsize=9)
    ax.legend(fontsize=7, loc='lower left'); plt.tight_layout(); plt.savefig(os.path.join(OUT, f'fig_D4_band_{c}.png'), dpi=110); plt.close()
    print(c, 'ref line %.1f km' % lref_km, L3[-1])
for fn, R in (('D4_line_displacement.csv', L1), ('D4_reach_displacement.csv', L2), ('D4_band_summary.csv', L3), ('D4_band_reaches.csv', L4)):
    keys = list(dict.fromkeys(kk for r in R for kk in r))
    with open(os.path.join(OUT, fn), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(R)
# summary figure: median and P90 |displacement| per single factor, per case
labs = [l for l, _ in FACT] + ['FRMP EL11 T50']
fig, ax = plt.subplots(figsize=(12, 4.5)); w = 0.2
for i, c in enumerate(CASES):
    med = []; p90 = []
    for lab, key in FACT + [('FRMP EL11 T50', None)]:
        nm = 'FRMP EL11 T50' if key is None else skey(key)
        r = [x for x in L1 if x['case'] == c and x['scenario'] == nm]
        med.append(r[0]['median_abs_disp_m'] if r else np.nan); p90.append(r[0]['p90_abs_disp_m'] if r else np.nan)
    x = np.arange(len(labs)) + (i - 1.5) * w
    ax.bar(x, med, w, label=f"{c} {CASES[c]['name']}"); ax.scatter(x, p90, s=10, c='k', zorder=3)
ax.set_xticks(np.arange(len(labs))); ax.set_xticklabels(labs, rotation=20, ha='right', fontsize=8); ax.set_yscale('symlog', linthresh=10)
ax.set_ylabel('shift of the T50 flood line (m)\nbars = median, dots = 90th percentile'); ax.legend(fontsize=8); ax.grid(axis='y', alpha=0.3)
ax.set_title('D4 – how far the T50 flood line (C1, CN II) moves when one assumption changes', fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'fig_D4_displacement.png'), dpi=120)
for r in L1: print(r['case'], r['scenario'].ljust(22), r['median_abs_disp_m'], r['p90_abs_disp_m'], r['share_gt25m'], r['mean_band_width_m'], r['bld_newly_in'], r['bld_newly_out'])
