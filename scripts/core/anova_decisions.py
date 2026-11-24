"""Variance decomposition of the decision metrics over the factorial of unknowns (Paper 1, section 3.5).

Two-way ANOVA without replication (SS of factor A, factor B, and A x B interaction = residual) for
four designs per case, all on the existing HEC-RAS max-depth rasters:
  T   : channel {C0, C1, C+} x return period {T50, T100}         (CN II, base maintenance)
  AMC : channel {C0, C1, C+} x antecedent moisture {II, III}     (T50, base)
  H   : channel {C0, C1, C+} x hydrology {T50, T100, T50 CN III} (base)
  M   : channel {C1, C+}     x maintenance {base, VegOver, VegClear, Block50, Sed05}  (T50 CN II)
Responses
  scalar : A010 (ha, depth >= 0.10 m in the domain), bld010, bld050 (buildings with max depth >= 0.10 / 0.50 m)
  spatial: wet indicator of every pixel (where the line falls, D4) and of every building (which buildings are in, D1/D4);
           SS summed over pixels / buildings -> share of the 'disagreement' between runs explained by each factor.
Rain thresholds (D3): channel x AMC x duration on ln(P*) for the 50-building threshold (D3_rain_thresholds.csv),
           three-way ANOVA without replication (main effects + 2-way interactions; 3-way = residual).

usage: P1_ROOT=<02_CASES> python anova_decisions.py   (writes SER-P1_paper_decisions/ANOVA_decision_metrics/)
"""
import os, sys, csv, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p1_common import ROOT, P1, CASES, Case, WET, DEEP

OUT = os.path.join(P1, 'ANOVA_decision_metrics'); os.makedirs(OUT, exist_ok=True)

DESIGNS = {
    'T':   ('channel', ['C0', 'C1', 'C+'], 'return period', ['T50', 'T100'],
            lambda a, b: (a, b, 'II', 'base')),
    'AMC': ('channel', ['C0', 'C1', 'C+'], 'AMC', ['II', 'III'],
            lambda a, b: (a, 'T50', b, 'base')),
    'H':   ('channel', ['C0', 'C1', 'C+'], 'hydrology', ['T50', 'T100', 'T50 CN III'],
            lambda a, b: (a, 'T50', 'III', 'base') if b == 'T50 CN III' else (a, b, 'II', 'base')),
    'M':   ('channel', ['C1', 'C+'], 'maintenance', ['base', 'VegOver', 'VegClear', 'Block50', 'Sed05'],
            lambda a, b: (a, 'T50', 'II', b)),
}

def ss2(Y):
    """Y: (nA, nB, ...) -> SS_A, SS_B, SS_AB, SS_T summed over trailing axes"""
    Y = Y.astype('float64'); m = Y.mean(axis=(0, 1), keepdims=True)
    mA = Y.mean(axis=1, keepdims=True); mB = Y.mean(axis=0, keepdims=True)
    nA, nB = Y.shape[:2]
    SSA = (nB * (mA - m) ** 2).sum(); SSB = (nA * (mB - m) ** 2).sum()
    SST = ((Y - m) ** 2).sum(); return SSA, SSB, SST - SSA - SSB, SST

rows = []
for c in CASES:
    K = Case(c); b, lab = K.buildings(); cache = {}
    def get(key):
        if key not in cache:
            D = K.depth(key); bd = K.bld_depth(D)
            cache[key] = dict(wet=(D >= WET)[K.dom], bw=bd >= WET,
                              A010=float((D >= WET).sum() * K.pxa), bld010=int((bd >= WET).sum()), bld050=int((bd >= DEEP).sum()))
        return cache[key]
    for dn, (fa, la, fb, lb, keyf) in DESIGNS.items():
        runs = [[get(keyf(a, x)) for x in lb] for a in la]
        for resp in ['A010', 'bld010', 'bld050', 'pixels', 'buildings']:
            if resp == 'pixels':
                Y = np.array([[r['wet'] for r in rr] for rr in runs])
            elif resp == 'buildings':
                Y = np.array([[r['bw'] for r in rr] for rr in runs])
            else:
                Y = np.array([[r[resp] for r in rr] for rr in runs], dtype=float)
            SSA, SSB, SSI, SST = ss2(Y)
            vals = Y.reshape(len(la), len(lb), -1).mean(axis=2) if Y.ndim > 2 else Y
            rows.append(dict(case=c, design=dn, factor_A=fa, factor_B=fb, response=resp,
                             share_A=round(SSA / SST, 3) if SST > 0 else np.nan,
                             share_B=round(SSB / SST, 3) if SST > 0 else np.nan,
                             share_AxB=round(SSI / SST, 3) if SST > 0 else np.nan,
                             range_over_runs=(round(float(vals.min()), 1), round(float(vals.max()), 1)) if Y.ndim == 2 else '',
                             n_units=int(np.prod(Y.shape[2:])) if Y.ndim > 2 else 1))
        print(c, dn, 'done', flush=True)

with open(os.path.join(OUT, 'anova_runs.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

# ---- rain thresholds (D3): ln P* for 50 buildings, channel x AMC x duration
rt = list(csv.DictReader(open(os.path.join(P1, 'D3_thresholds', 'D3_rain_thresholds.csv'))))
DUR = ['P_1h', 'P_2h', 'P_3h', 'P_6h', 'P_12h', 'P_24h']
rr = []
for c in CASES:
    sub = [r for r in rt if r['case'] == c and r['threshold'] == '50 buildings' and r['P_1h'] not in ('', None)]
    chans = sorted({r['channel'] for r in sub}, key=['C0', 'C1', 'C+'].index)
    if len(chans) < 2 or len({r['P_1h'] for r in sub}) == 1 and len(chans) == 1:
        continue
    Y = np.array([[[np.log(float(next(r for r in sub if r['channel'] == ch and r['AMC'] == am)[d])) for d in DUR]
                   for am in ['AMCII', 'AMCIII']] for ch in chans])           # (ch, AMC, D)
    m = Y.mean(); SST = ((Y - m) ** 2).sum()
    out = dict(case=c, channels='/'.join(chans))
    names = ['channel', 'AMC', 'duration']
    tot = 0
    for i, n in enumerate(names):
        ax = tuple(j for j in range(3) if j != i)
        mi = Y.mean(axis=ax, keepdims=True); s = ((mi - m) ** 2).sum() * Y.size / Y.shape[i]
        out['share_' + n] = round(s / SST, 3); tot += s
    for i, j in itertools.combinations(range(3), 2):
        k = 3 - i - j
        mij = Y.mean(axis=k, keepdims=True); mi = Y.mean(axis=(j, k), keepdims=True); mj = Y.mean(axis=(i, k), keepdims=True)
        s = ((mij - mi - mj + m) ** 2).sum() * Y.shape[k]
        out[f'share_{names[i]}x{names[j]}'] = round(s / SST, 3); tot += s
    out['share_residual'] = round(1 - tot / SST, 3)
    rr.append(out)
if rr:
    with open(os.path.join(OUT, 'anova_rain_thresholds.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rr[0])); w.writeheader(); w.writerows(rr)
print('written', OUT)
