"""Paper 1, Fig. 4 – D2 is upkeep worth it: change in flooded area (%) and in buildings reached (>= 0.10 m) for each maintenance state vs
the same channel unmaintained (T50, CN II), per torrent and channel state.
Input: ../tables/raw/scen_D2_summary.csv.  usage: python fig04_d2_maintenance.py"""
import os, csv
import numpy as np
from p1_style import plt, RAW, CASES, NAMES, CH_COL, save
R = list(csv.DictReader(open(os.path.join(RAW, 'scen_D2_summary.csv'))))
SC = [('VegOver', 'overgrown\n(n 0.08)'), ('VegClear', 'cleared\n(n 0.03)'), ('Block50', 'crossings\n50% blocked'), ('Sed05', 'bed\n+0.5 m')]
fig, axs = plt.subplots(2, 4, figsize=(13, 6.2), sharex=True)
for j, c in enumerate(CASES):
    for i, (resp, lab) in enumerate([('dA_pct', 'Δ flooded area (%)'), ('dbld', 'Δ buildings ≥ 0.10 m')]):
        ax = axs[i, j]
        for k, ch in enumerate(['C1', 'Cp']):
            vals = []
            for s, _ in SC:
                r = next(x for x in R if x['case'] == c and x['channel'] == ch and x['scenario'] == s)
                vals.append(float(r['dA_pct']) if resp == 'dA_pct' else int(r['bld010']) - int(r['bld010_base']))
            x = np.arange(4) + (k - .5) * .38
            ax.bar(x, vals, width=.36, color=CH_COL['C+' if ch == 'Cp' else ch], label=('C1 assumed channel' if ch == 'C1' else 'C+ larger channel'))
            for xx, v in zip(x, vals):
                if abs(v) >= (40 if resp == 'dA_pct' else 0.25 * max(1, max(abs(u) for u in vals))):
                    ax.text(xx, v, f'{v:+.0f}', ha='center', va='bottom' if v > 0 else 'top', fontsize=7)
        ax.axhline(0, color='#333', lw=.8)
        ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l in SC], fontsize=7.5)
        if j == 0: ax.set_ylabel(lab)
        if i == 0:
            b = [r for r in R if r['case'] == c and r['channel'] == 'C1'][0]
            ax.set_title(f'{c} {NAMES[c]}\nbase C1: {float(b["A_base_ha"]):.0f} ha, {b["bld010_base"]} buildings', loc='left', fontsize=9)
axs[0, 0].legend(loc='upper right', frameon=False)
lo = min(ax.get_ylim()[0] for ax in axs[0]); hi = max(ax.get_ylim()[1] for ax in axs[0])
for ax in axs[0]: ax.set_ylim(lo, hi)
for ax in axs[1]:
    a, b = ax.get_ylim(); ax.set_ylim(min(a, -60), max(b, 60))   # at least +-60 buildings, so small changes look small          # area: common scale; buildings: own scale per torrent (base counts differ 59-2156)
fig.tight_layout(); save(fig, 'Fig04_D2_maintenance.png')
