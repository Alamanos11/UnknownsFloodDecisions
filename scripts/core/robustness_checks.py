"""Paper 1 – robustness checks prompted by the internal review (24/11/2026).
(1) band of plausible T50 lines restricted to the states a channel survey can resolve (CN II: channel x maintenance, 11 states)
    and to the channel alone (C0/C1/C+ at CN II, 3 states);
(2) variance shares of the T and AMC designs with two channel levels only ({C1, C+} and {C0, C1});
usage: P1_ROOT=<02_CASES> python robustness_checks.py <out_dir>   (writes D4_band_ensembles.csv = robust_band.csv and anova_channel_levels.csv)"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p1_common import CASES, MAINT, Case, WET, DEEP
import scipy.ndimage as ndi
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
def clean(E, k):
    E = ndi.binary_opening(E, iterations=1); E = ndi.binary_fill_holes(E)
    lab, n = ndi.label(E)
    if n:
        sz = ndi.sum(E, lab, range(1, n + 1)) * k.px * k.px; E = np.isin(lab, 1 + np.where(sz >= 400)[0])
    r = max(1, int(round(4 / k.px))); st = ndi.iterate_structure(ndi.generate_binary_structure(2, 1), r)
    return ndi.binary_closing(np.pad(E, r), structure=st)[r:-r, r:-r] & k.dom
def ss2(Y):
    Y = Y.astype('float64'); m = Y.mean(axis=(0, 1), keepdims=True)
    mA = Y.mean(axis=1, keepdims=True); mB = Y.mean(axis=0, keepdims=True); nA, nB = Y.shape[:2]
    SSA = (nB * (mA - m) ** 2).sum(); SSB = (nA * (mB - m) ** 2).sum(); SST = ((Y - m) ** 2).sum()
    return SSA / SST, SSB / SST, 1 - (SSA + SSB) / SST
ENS = {'14 states (paper)': [(ch, 'T50', cn, 'base') for ch in ('C0', 'C1', 'C+') for cn in ('II', 'III')] + [(ch, 'T50', 'II', m) for ch in ('C1', 'C+') for m in MAINT],
       '11 states, CN II (channel x maintenance)': [(ch, 'T50', 'II', 'base') for ch in ('C0', 'C1', 'C+')] + [(ch, 'T50', 'II', m) for ch in ('C1', 'C+') for m in MAINT],
       '3 states, CN II (channel only)': [(ch, 'T50', 'II', 'base') for ch in ('C0', 'C1', 'C+')],
       '2 states, CN II (C1, C+)': [(ch, 'T50', 'II', 'base') for ch in ('C1', 'C+')]}
R1, R2 = [], []
for c in CASES:
    k = Case(c); b, _ = k.buildings(); cache = {}
    def E(key):
        if key not in cache:
            D = k.depth(key); cache[key] = dict(E=clean(D >= WET, k), wet=(D >= WET)[k.dom], bw=k.bld_depth(D) >= WET, bd=k.bld_depth(D) >= DEEP)
        return cache[key]
    for name, ens in ENS.items():
        ens = [x for x in ens if x in k.cf['plans']]; N = np.zeros(k.sh, 'uint8'); BN = np.zeros(len(b), 'uint8')
        for key in ens:
            e = E(key)['E']; N += e; BN += (k.bld_depth(e.astype('float32')) > 0)
        n = len(ens)
        R1.append(dict(case=c, ensemble=name, states=n, band_ha=round((((N > 0) & (N < n)) & k.dom).sum() * k.pxa, 1),
                       band_buildings=int(((BN > 0) & (BN < n)).sum()), always_ha=round(((N == n) & k.dom).sum() * k.pxa, 1)))
    for chans in (['C1', 'C+'], ['C0', 'C1'], ['C0', 'C1', 'C+']):
        for dn, lev, keyf in (('T', ['T50', 'T100'], lambda a, x: (a, x, 'II', 'base')), ('AMC', ['II', 'III'], lambda a, x: (a, 'T50', x, 'base'))):
            runs = [[E(keyf(a, x)) for x in lev] for a in chans]
            for resp, f in (('pixels', 'wet'), ('buildings', 'bw'), ('bld050', 'bd')):
                Y = np.array([[r[f] for r in rr] for rr in runs])
                if resp == 'bld050': Y = Y.sum(axis=2)
                sa, sb, si = ss2(Y)
                R2.append(dict(case=c, design=dn, channels='/'.join(chans), response=resp, share_channel=round(sa, 3), share_other=round(sb, 3), share_interaction=round(si, 3)))
    print(c, 'done', flush=True)
for nm, R in (('robust_band.csv', R1), ('robust_anova_channel_levels.csv', R2)):
    with open(os.path.join(OUT, nm), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(R[0])); w.writeheader(); w.writerows(R)
print('written', OUT)
