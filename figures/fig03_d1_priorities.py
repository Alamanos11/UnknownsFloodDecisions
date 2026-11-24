"""Paper 1, Fig. 3 – D1 which stream / reach first: (a) stream order over all mixed-state combinations; (b) probability that a reach is
in the top 10 (by buildings reached) over the mixed states; (c) Spearman rho of the reach ranking vs the reference (C1 T50) per scenario.
Input: ../tables/raw/D1_stream_order.csv, D1_mixed_states.csv, D1_uniform_scenarios.csv.  usage: python fig03_d1_priorities.py"""
import os, csv
import numpy as np
from p1_style import plt, RAW, CASE_COL, save
so = list(csv.DictReader(open(os.path.join(RAW, 'D1_stream_order.csv'))))
ms = list(csv.DictReader(open(os.path.join(RAW, 'D1_mixed_states.csv'))))
us = list(csv.DictReader(open(os.path.join(RAW, 'D1_uniform_scenarios.csv'))))
fig = plt.figure(figsize=(12, 5.2)); gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1.5, 1.1], height_ratios=[1, 1.6])
# (a)
ax = fig.add_subplot(gs[0, 0])
orders = sorted({r['stream_order'] for r in so}, key=lambda o: -max(float(r['share']) for r in so if r['stream_order'] == o))
shade = ['#2a78d6', '#9cc3ee', '#d7e7f8']
for y, cn in enumerate(['II', 'III']):
    left = 0
    for k, o in enumerate(orders):
        r = next((x for x in so if x['CN'] == cn and x['stream_order'] == o), None)
        if r is None: continue
        v = float(r['share']); ax.barh(y, v, left=left, color=shade[k], edgecolor='white', linewidth=1.5, height=.6)
        if v > .08: ax.text(left + v / 2, y, f'{v * 100:.0f}%', ha='center', va='center', fontsize=8, color='white' if k == 0 else '#222')
        left += v
ax.set_yticks([0, 1]); ax.set_yticklabels(['CN II\n(14 641 comb.)', 'CN III\n(81 comb.)']); ax.set_xlim(0, 1); ax.invert_yaxis()
ax.set_xticks([0, .5, 1]); ax.set_xticklabels(['0', '50', '100%']); ax.grid(False)
for k, o in enumerate(orders): ax.bar(0, 0, color=shade[k], label=o.replace(' > ', ' > '))
ax.legend(loc='upper center', bbox_to_anchor=(.5, -.25), fontsize=7.2, frameon=False)
ax.set_title('(a) Order of the four streams\n(each stream in any of its states)', loc='left')
# (b)
ax = fig.add_subplot(gs[:, 1])
rr = sorted([r for r in ms if r['CN'] == 'II' and float(r['P_top10']) > 0], key=lambda r: -float(r['P_top10']))[:20]
y = np.arange(len(rr))
ax.barh(y, [float(r['P_top10']) for r in rr], color=[CASE_COL[r['reach'][:6]] for r in rr], height=.7)
ax.set_yticks(y); ax.set_yticklabels([f"{r['reach']} (ch {r['ch_from']}–{r['ch_to']} m)" for r in rr], fontsize=7.5); ax.invert_yaxis()
ax.axvline(.9, color='#555', ls=':', lw=1); ax.text(.9, -.9, '0.9', fontsize=7.5, color='#555', ha='center')
ax.set_xlim(0, 1); ax.set_xlabel('P(reach in top 10), CN II mixed states')
for c, col in CASE_COL.items(): ax.barh(0, 0, color=col, label=c)
ax.legend(loc='lower right', frameon=False)
ax.set_title('(b) Which reaches stay in the top 10 (71 reaches)', loc='left')
# (c)
ax = fig.add_subplot(gs[:, 2])
sc = [r for r in us if r['scenario'] != 'C1 T50']
v = np.array([float(r['rho_all_reaches']) for r in sc]); grp = ['maintenance' if any(m in r['scenario'] for m in ['Veg', 'Block', 'Sed']) else 'channel / hydrology' for r in sc]
for g, col, x0 in [('channel / hydrology', '#2a78d6', 0), ('maintenance', '#eb6834', 1)]:
    vv = v[[i for i, x in enumerate(grp) if x == g]]
    ax.scatter(x0 + np.linspace(-.18, .18, len(vv)), vv, color=col, s=26, zorder=3)
    ax.hlines(np.median(vv), x0 - .25, x0 + .25, color='#222', lw=1.4)
ax.set_xticks([0, 1]); ax.set_xticklabels(['channel /\nhydrology (9)', 'maintenance\n(8)']); ax.set_ylim(.6, 1.0)
ax.set_ylabel('Spearman ρ vs reference (C1 T50)'); ax.set_xlim(-.6, 1.6)
ax.set_title('(c) Reach ranking per scenario\n(Spearman ρ, 71 reaches)', loc='left')
fig.tight_layout(); save(fig, 'Fig03_D1_priorities.png')
