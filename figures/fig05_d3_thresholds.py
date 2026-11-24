"""Paper 1, Fig. 5 – D3 when damage starts: outlet discharge at which 10 / 50 / 100 buildings are reached (>= 0.10 m), per torrent and
channel state, as a fraction of the T50 peak (stepped-inflow HEC-RAS runs). Open markers with arrows = below the first step (upper bound).
Input: ../tables/raw/D3_thresholds_summary.csv.  usage: python fig05_d3_thresholds.py"""
import os, csv
import numpy as np
from p1_style import plt, RAW, CASES, NAMES, CH_COL, CH_LAB, save
R = list(csv.DictReader(open(os.path.join(RAW, 'D3_thresholds_summary.csv'))))
fig, ax = plt.subplots(figsize=(8.5, 4.6))
MK = {'Q_10bld': ('v', '10'), 'Q_50bld': ('o', '50'), 'Q_100bld': ('^', '100')}
ylab = []
for i, c in enumerate(CASES):
    for k, ch in enumerate(['C0', 'C1', 'Cp']):
        r = next(x for x in R if x['case'] == c and x['channel'] == ch)
        y = i * 4 + k; q50 = float(r['Q_T50']); chk = 'C+' if ch == 'Cp' else ch
        xs = []
        for col, (m, _) in MK.items():
            v = r[col]; up = v.startswith('<'); q = float(v.lstrip('<')) / q50; xs.append(q)
            ax.scatter(q, y, marker=m, s=38, facecolor='white' if up else CH_COL[chk], edgecolor=CH_COL[chk], zorder=3, lw=1.2)
            if up: ax.annotate('', xy=(q * .8, y), xytext=(q, y), arrowprops=dict(arrowstyle='->', color=CH_COL[chk], lw=1))
        ax.plot([min(xs), max(xs)], [y, y], color=CH_COL[chk], lw=1, zorder=2)
        ylab.append((y, f'{c} {chk}'))
    ax.text(3.05, i * 4 + 1, f'T50 = {float(r["Q_T50"]):.0f} m³/s', fontsize=7.5, va='center', color='#555')
ax.axvline(1, color='#333', lw=.9, ls='--'); ax.text(1.02, -1.2, 'T50 peak', fontsize=7.5)
ax.set_xscale('log'); ax.set_xlim(.08, 3); ax.set_xticks([.1, .2, .5, 1, 2]); ax.set_xticklabels(['0.1', '0.2', '0.5', '1', '2'])
ax.set_yticks([y for y, _ in ylab]); ax.set_yticklabels([l for _, l in ylab], fontsize=7.8); ax.invert_yaxis()
ax.set_xlabel('outlet discharge at which N buildings are reached / T50 peak discharge')
for col, (m, n) in MK.items(): ax.scatter([], [], marker=m, color='#555', label=f'{n} buildings')
ax.scatter([], [], marker='o', facecolor='white', edgecolor='#555', label='below first step (upper bound)')
ax.legend(loc='lower right', fontsize=7.5, frameon=True)
fig.tight_layout(); save(fig, 'Fig05_D3_thresholds.png')
