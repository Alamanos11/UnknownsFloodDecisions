"""Shared style for the Paper 1 figures (journal: Journal of Hydrology; 300 dpi PNG; colour-vision-safe palette)."""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, '..', 'tables', 'raw')          # raw result CSVs (copied from 02_CASES)
CASES = ['SER-01', 'SER-02', 'SER-03', 'SER-04']
NAMES = {'SER-01': 'Agioi Anargyroi', 'SER-02': 'Agia Varvara', 'SER-03': 'Agios Ioannis', 'SER-04': 'Krousovitis'}
CASE_COL = {'SER-01': '#6a3d9a', 'SER-02': '#b15928', 'SER-03': '#008080', 'SER-04': '#d6457f'}   # cases (purple, brown, teal, magenta)
CH_COL = {'C0': '#e07b39', 'C1': '#2a78d6', 'C+': '#1baf7a'}                                       # channel states
CH_LAB = {'C0': 'C0 DEM only', 'C1': 'C1 assumed channel', 'C+': 'C+ larger channel'}
plt.rcParams.update({'font.size': 8.5, 'axes.titlesize': 9.5, 'axes.labelsize': 8.5, 'legend.fontsize': 8,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 300, 'figure.dpi': 100,
                     'axes.grid': True, 'grid.color': '#e6e6e3', 'grid.linewidth': .6, 'axes.axisbelow': True})
def save(fig, name):
    out = os.path.join(HERE, name); fig.savefig(out, dpi=300, bbox_inches='tight'); print('saved', out)
