"""Paper 1, Fig. 2 – study design: unknowns -> HEC-HMS / HEC-RAS 2D runs -> decision metrics -> decision answers.
usage: python fig02_design.py"""
import sys
from p1_style import plt, save
from matplotlib.patches import FancyBboxPatch
fig, ax = plt.subplots(figsize=(13, 6.6)); ax.set_xlim(0, 13); ax.set_ylim(0, 6.6); ax.axis('off')
def box(x, y, w, h, title, lines, fc, ec):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02,rounding_size=0.12', fc=fc, ec=ec, lw=1.3))
    ax.text(x + .15, y + h - .18, title, fontsize=10, fontweight='bold', va='top')
    ax.text(x + .15, y + h - .55, '\n'.join(lines), fontsize=8.6, va='top', linespacing=1.45)
def arrow(x0, y0, x1, y1):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle='-|>', color='#444', lw=1.4))
ax.text(.1, 6.45, 'Unknowns (no survey, no gauge)', fontsize=9, color='#555')
box(.1, 3.7, 2.9, 2.6, 'Epistemic unknowns', ['Channel: C0 DEM only / C1 assumed', '   / C+ larger', 'Soil moisture: AMC II / III',
    'Rainfall product (events):', '   gauge vs reanalysis'], '#e8f1fb', '#2a78d6')
box(.1, .6, 2.9, 2.75, 'Choices and levers', ['Design standard: T50 / T100', '   (+ T1000)', 'Maintenance state:', '   overgrown, cleared,',
    '   half-blocked crossings,', '   bed +0.5 m'], '#fdeee6', '#eb6834')
ax.text(3.45, 6.45, 'Models (4 torrents)', fontsize=9, color='#555')
box(3.45, .6, 2.8, 5.7, 'HEC-HMS + HEC-RAS 2D', ['National IDF, SCS-CN,', 'alternating-block storms', '',
    '46 core and event runs', '44 maintenance and', 'stepped-inflow runs', '', 'one mesh per torrent', '(6–55 km², 61–122 k cells)',
    'closed-edge check', 'on every run', '', 'Buildings: OSM + MS', 'footprints'], '#f3f3f1', '#777')
ax.text(6.7, 6.45, 'Decision metrics', fontsize=9, color='#555')
box(6.7, 4.75, 2.9, 1.55, 'D1 Which stream first', ['rank of streams and 71 reaches', 'by buildings reached; ρ, P(top-10)'], '#ffffff', '#2a78d6')
box(6.7, 3.15, 2.9, 1.45, 'D2 Is upkeep worth it', ['Δ area, Δ buildings vs the same', 'channel unmaintained'], '#ffffff', '#2a78d6')
box(6.7, 1.6, 2.9, 1.4, 'D3 When to warn', ['outlet Q for 10/50/100 buildings', '(stepped inflow to Q1000)'], '#ffffff', '#2a78d6')
box(6.7, .05, 2.9, 1.4, 'D4 Where the line falls', ['line shift (median, P90);', 'band over 14 states + buildings'], '#ffffff', '#2a78d6')
ax.text(10.05, 6.45, 'Analysis and answers', fontsize=9, color='#555')
box(10.05, 3.15, 2.85, 3.15, 'Which unknown matters?', ['two-way ANOVA per torrent:', 'channel × T, channel × AMC,', 'channel × maintenance;',
    'pixel, building and', 'scalar responses'], '#e9f7f1', '#1baf7a')
box(10.05, .05, 2.85, 2.95, 'Benchmarks (no calibration)', ['official FRMP maps', '(extent, depth)', 'Sentinel-1 change detection',
    'Jan 2021 events:', 'gauge vs reanalysis rain'], '#f3f3f1', '#777')
arrow(3.0, 5.0, 3.45, 5.0); arrow(3.0, 2.0, 3.45, 2.0)
for y in (5.5, 3.85, 2.3, .75): arrow(6.25, y, 6.7, y)
arrow(9.6, 4.6, 10.05, 4.6); arrow(9.6, 1.5, 10.05, 1.5)
save(fig, 'Fig02_design.png')
