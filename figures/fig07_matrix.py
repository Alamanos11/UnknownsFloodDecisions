"""Paper 1, Fig. 7 – decision x unknown matrix: share (%) of the between-run variance of each decision metric explained by
the channel unknown and by the second unknown in three two-way designs (anova_runs.csv from anova_decisions.py).
usage: python fig07_matrix.py [anova_runs.csv]   (default ../tables/raw/anova_runs.csv; from anova_decisions.py)"""
import sys, csv
import numpy as np
import os
from p1_style import plt, RAW, save
plt.rcParams['axes.grid'] = False
R = list(csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(RAW, 'anova_runs.csv'))))
CASES = ['SER-01', 'SER-02', 'SER-03', 'SER-04']
RESP = [('pixels', 'Where the line falls\n(D4; wet/dry per pixel)'), ('buildings', 'Which buildings are in\n(D1, D4; per building)'),
        ('bld050', 'Damaging flooding\n(buildings ≥ 0.5 m)'), ('A010', 'Flooded area\n(≥ 0.10 m)')]
DES = [('T', 'Design standard\nT50 vs T100'), ('AMC', 'Soil moisture\nAMC II vs III'), ('M', 'Maintenance state\n(5 states, C1/C+)')]
M = np.zeros((len(RESP) * len(CASES), len(DES) * 3))
for i, (rk, _) in enumerate(RESP):
    for j, c in enumerate(CASES):
        for k, (dn, _) in enumerate(DES):
            r = next(x for x in R if x['case'] == c and x['design'] == dn and x['response'] == rk)
            M[i * 4 + j, k * 3:(k + 1) * 3] = [float(r['share_A']), float(r['share_B']), float(r['share_AxB'])]
fig, ax = plt.subplots(figsize=(10.5, 9.2))
im = ax.imshow(M * 100, cmap='Blues', vmin=0, vmax=100, aspect='auto')
for (y, x), v in np.ndenumerate(M):
    ax.text(x, y, f'{v * 100:.0f}', ha='center', va='center', fontsize=8.5, color='white' if v > .6 else '#222')
ax.set_xticks(range(9)); ax.set_xticklabels(['channel', 'T50/T100', 'interaction', 'channel', 'AMC', 'interaction',
                                             'channel', 'maintenance', 'interaction'], rotation=35, ha='right', fontsize=8.5)
ax.set_yticks(range(16)); ax.set_yticklabels(CASES * 4, fontsize=8.5)
for k, (_, t) in enumerate(DES):
    ax.text(k * 3 + 1, -1.0, t, ha='center', va='bottom', fontsize=9.5, fontweight='bold')
    if k: ax.axvline(k * 3 - .5, color='white', lw=4)
for i, (_, t) in enumerate(RESP):
    ax.text(-1.6, i * 4 + 1.5, t, ha='right', va='center', fontsize=9)
    if i: ax.axhline(i * 4 - .5, color='white', lw=4)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=.03, pad=.02); cb.set_label('share of the variance between runs (%)', fontsize=8.5)
fig.subplots_adjust(left=.27, right=.93, top=.88, bottom=.1)
save(fig, 'Fig07_decision_matrix.png')
