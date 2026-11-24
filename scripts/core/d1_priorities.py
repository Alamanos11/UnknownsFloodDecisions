"""D1 - which streams / reaches to study first: is the priority order robust to the unknowns?

Exposure = buildings reached by >= 0.10 m (and >= 0.50 m) of max depth anywhere on their footprint, per reach.
A building belongs to the reach of the nearest centreline station (reaches between crossings, <= 500 m, as in D3).
Reference ranking: C1, T50, CN II, channel as assumed. Tests:
  1. uniform scenarios (the same state applied to all four streams): Spearman rho of the reach exposures
     against the reference, overlap of the top-10 reaches, order of the four streams;
  2. mixed states: every combination of an independent channel / maintenance state per stream
     (CN II: C0, C1, C+, C1 x 4 maintenance, C+ x 4 maintenance = 11 states -> 11^4 = 14 641 combinations;
      CN III: C0, C1, C+ -> 3^4 = 81 combinations; soil moisture is regional, so CN is common to all streams):
     probability that each reach is in the regional top 10, and that each stream ranks first.
Writes SER-P1_paper_decisions/data/reaches.gpkg (reach segments, for maps) and SER-P1_paper_decisions/D1_priorities/: D1_reach_exposure.csv, D1_uniform_scenarios.csv, D1_mixed_states.csv,
D1_stream_order.csv, fig_D1_reaches.png.
usage: P1_ROOT=<02_CASES> python3 d1_priorities.py"""
import os, itertools, csv, numpy as np
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from p1_common import CASES, REF, MAINT, Case, skey, P1, WET, DEEP

OUT = os.path.join(P1, 'D1_priorities'); os.makedirs(OUT, exist_ok=True)
rows, R = [], {}          # R[case][key] = (bld010 per reach, bld050 per reach)
reaches = []; RL = []
for c in CASES:
    k = Case(c); b, _ = k.buildings(); nr = len(k.edges) - 1; rt = k.reach_table(); reaches += rt; RL.append(k.reach_lines())
    nb_reach = np.bincount(b['reach_i'], minlength=nr)
    R[c] = {}
    for key in k.cf['plans']:
        D = k.depth(key); bd = k.bld_depth(D)
        e1 = np.bincount(b['reach_i'][bd >= WET], minlength=nr); e5 = np.bincount(b['reach_i'][bd >= DEEP], minlength=nr)
        A = np.bincount(k.reach[(D >= WET) & k.dom], minlength=nr) * k.pxa
        R[c][key] = (e1, e5)
        for j in range(nr):
            rows.append(dict(**rt[j], buildings_in_reach=int(nb_reach[j]), scenario=skey(key), channel=key[0], T=key[1], CN=key[2],
                             maintenance=key[3], bld010=int(e1[j]), bld050=int(e5[j]), A_wet_ha=round(float(A[j]), 2)))
    print(c, len(b), 'buildings,', nr, 'reaches,', len(k.cf['plans']), 'scenarios')
import pandas as pd, geopandas as gpd
gpd.GeoDataFrame(pd.concat(RL, ignore_index=True), crs=RL[0].crs).to_file(os.path.join(P1, 'data', 'reaches.gpkg'))
with open(os.path.join(OUT, 'D1_reach_exposure.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

names = [r['reach'] for r in reaches]; cases = list(CASES)
def vec(states, i=0):   # concatenated reach exposure for one state per case
    return np.concatenate([R[c][s][i] for c, s in zip(cases, states)])
def top(v, n=10): return set(np.argsort(-v, kind='stable')[:n][v[np.argsort(-v, kind='stable')[:n]] > 0])
def order(states, i=0):
    t = [R[c][s][i].sum() for c, s in zip(cases, states)]; return [cases[j] for j in np.argsort(-np.array(t), kind='stable')], t
ref = vec([REF] * 4); ref5 = vec([REF] * 4, 1); top_ref = top(ref)

# 1. uniform scenarios
common = [k for k in CASES['SER-01']['plans'] if all(k in CASES[c]['plans'] for c in cases)]
U = []
for key in common:
    v = vec([key] * 4); v5 = vec([key] * 4, 1); o, t = order([key] * 4)
    win = {c: (spearmanr(R[c][key][0], R[c][REF][0])[0] if R[c][REF][0].std() > 0 and R[c][key][0].std() > 0 else np.nan) for c in cases}
    U.append(dict(scenario=skey(key), rho_all_reaches=round(spearmanr(v, ref)[0], 3), rho_050=round(spearmanr(v5, ref5)[0], 3),
                  top10_overlap=len(top(v) & top_ref), stream_order=' > '.join(o), **{f'bld010_{c}': int(x) for c, x in zip(cases, t)},
                  **{f'rho_within_{c}': round(win[c], 3) for c in cases},
                  **{f'top_reach_{c}': names[sum(len(R[cc][REF][0]) for cc in cases[:i]) + int(np.argmax(R[c][key][0]))] for i, c in enumerate(cases)}))
with open(os.path.join(OUT, 'D1_uniform_scenarios.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(U[0])); w.writeheader(); w.writerows(U)

# 2. mixed states
ST = {'II': [('C0', 'T50', 'II', 'base'), ('C1', 'T50', 'II', 'base'), ('C+', 'T50', 'II', 'base')] +
            [(ch, 'T50', 'II', m) for ch in ('C1', 'C+') for m in MAINT],
      'III': [('C0', 'T50', 'III', 'base'), ('C1', 'T50', 'III', 'base'), ('C+', 'T50', 'III', 'base')]}
M, SO = [], []
for cn, st in ST.items():
    cnt10 = np.zeros(len(names)); first = {c: 0 for c in cases}; orders = {}; rhos = []; n = 0
    for combo in itertools.product(st, repeat=4):
        v = vec(combo); n += 1
        for j in top(v): cnt10[j] += 1
        o, _ = order(combo); first[o[0]] += 1; orders[' > '.join(o)] = orders.get(' > '.join(o), 0) + 1
        rhos.append(spearmanr(v, ref)[0])
    for j, nm in enumerate(names):
        M.append(dict(CN=cn, combinations=n, reach=nm, ch_from=reaches[j]['ch_from'], ch_to=reaches[j]['ch_to'], bld010_ref=int(ref[j]),
                      in_ref_top10=j in top_ref, P_top10=round(cnt10[j] / n, 3)))
    for o, m in sorted(orders.items(), key=lambda x: -x[1]):
        SO.append(dict(CN=cn, combinations=n, stream_order=o, share=round(m / n, 4)))
    print(f'CN {cn}: {n} combinations, rho vs reference median {np.median(rhos):.2f} (P5 {np.percentile(rhos, 5):.2f}),',
          'P(stream first):', {c: round(first[c] / n, 3) for c in cases})
with open(os.path.join(OUT, 'D1_mixed_states.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(M[0])); w.writeheader(); w.writerows(M)
with open(os.path.join(OUT, 'D1_stream_order.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(SO[0])); w.writeheader(); w.writerows(SO)

# figure: reference top-15 reaches, range over the uniform T50 states, P(top 10) under mixed states (CN II)
T50 = [k for k in common if k[1] == 'T50']
mat = np.array([vec([k] * 4) for k in T50]); p10 = {m['reach']: m['P_top10'] for m in M if m['CN'] == 'II'}
idx = [j for j in np.argsort(-mat.max(0)) if mat[:, j].max() > 0][:20]
fig, ax = plt.subplots(figsize=(12, 5)); x = np.arange(len(idx))
ax.bar(x, ref[idx], color='#264653', label='reference C1 T50 CN II')
ax.vlines(x, mat[:, idx].min(0), mat[:, idx].max(0), color='#e76f51', lw=3, alpha=0.7, label=f'range over {len(T50)} T50 states (channel, CN, maintenance)')
for xi, j in zip(x, idx): ax.text(xi, mat[:, j].max() * 1.02 + 5, f"{p10[names[j]]:.2f}", ha='center', fontsize=7, color='#2a9d8f')
ax.set_xticks(x); ax.set_xticklabels([f"{names[j]}\n{reaches[j]['ch_from']}–{reaches[j]['ch_to']}" for j in idx], rotation=90, fontsize=7)
ax.set_ylabel('buildings reached (≥ 0.10 m)'); ax.legend(fontsize=8)
ax.set_title('D1 – reaches with the most buildings reached; green numbers = probability of being in the regional top 10 under mixed states (CN II)', fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'fig_D1_reaches.png'), dpi=120)
for u in U: print(u['scenario'], u['rho_all_reaches'], u['top10_overlap'], u['stream_order'])
