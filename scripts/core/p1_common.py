"""Shared configuration and helpers for the Paper 1 decision analyses (D1, D4).

All paths are relative to the 02_CASES folder, given by the environment variable P1_ROOT
(default: the folder two levels above this script, i.e. 02_CASES when the script sits in
SER-P1_paper_decisions/scripts/).

Inputs used (all already in 02_CASES):
  <case>/03_hecras/results/<plan>/*maxdepth.tif   max depth of every HEC-RAS plan (core runs and D2 runs)
  SER-P1_paper_decisions/data/model_perimeters.json  2D-area perimeter of each case (EPSG:2100)
  SER-P1_paper_decisions/data/centrelines/<case>_centreline.gpkg  stream centreline used to build the models
  _data_buildings/buildings_merged_serres.gpkg      OSM + Microsoft building footprints (see its README)
  _data_frmp/frmp_el11_extents/frmp_<T>_<box>.tif  FRMP EL11 flood extents (0/1, 5 m)
"""
import os, glob, json
import numpy as np, rasterio, geopandas as gpd
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, mapping

ROOT = os.environ.get('P1_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
P1 = os.path.join(ROOT, 'SER-P1_paper_decisions')
WET = 0.10          # flood-line / wet threshold (m), as for the delivered flood lines
DEEP = 0.50         # 'damaging' depth for buildings (m)

# Scenario key = (channel, T, CN, maintenance). Reference for D1 and D4 = ('C1', 'T50', 'II', 'base').
# Plan folders per case. C0 = DEM with crossings opened, no channel; C1 = assumed channel; C+ = larger channel.
def _plans(c0, c1, cp, extra=None):
    d = {}
    for ch, pl in (('C0', c0), ('C1', c1), ('C+', cp)):
        for k, v in pl.items():
            T, cn = k
            d[(ch, T, cn, 'base')] = v
    for k, v in (extra or {}).items():
        d[k] = v
    return d

MAINT = ['VegOver', 'VegClear', 'Block50', 'Sed05']   # n 0.08 / n 0.03 / half waterway blocked at crossings / bed +0.5 m

CASES = {
 'SER-01': dict(dir='SER-01_AgioiAnargyroi', name='Agioi Anargyroi',
   plans=_plans({('T50', 'II'): 'p07_T50_G0open', ('T100', 'II'): 'p08_T100_G0open', ('T50', 'III'): 'p09_T50CN3_G0open'},
                {('T50', 'II'): 'p01_T50_G2', ('T100', 'II'): 'p02_T100_G2', ('T1000', 'II'): 'p03_T1000_G2', ('T50', 'III'): 'p04_T50CN3_G2'},
                {('T50', 'II'): 'p10_T50_G1', ('T100', 'II'): 'p11_T100_G1', ('T50', 'III'): 'p12_T50CN3_G1'},
                {**{('C1', 'T50', 'II', m): f'p{13+i}_T50_C1_{m}' for i, m in enumerate(MAINT)},
                 **{('C+', 'T50', 'II', m): f'p{17+i}_T50_Cp_{m}' for i, m in enumerate(MAINT)}}),
   line='SER-01_centreline.gpkg', off=0.0, srange=(40, 5084), half=9.5, frmp='S123',
   crossings=[('Meg. Alexandrou bridge', 2055), ('road', 3030), ('road', 3255), ('railway', 3875)]),
 'SER-02': dict(dir='SER-02_AgiaVarvara', name='Agia Varvara',
   plans=_plans({('T50', 'II'): 'p05_T50_G0', ('T100', 'II'): 'p06_T100_G0', ('T50', 'III'): 'p07_T50CN3_G0'},
                {('T50', 'II'): 'p01_T50_G1', ('T100', 'II'): 'p02_T100_G1', ('T1000', 'II'): 'p03_T1000_G1', ('T50', 'III'): 'p04_T50CN3_G1'},
                {('T50', 'II'): 'p08_T50_G2', ('T100', 'II'): 'p09_T100_G2', ('T50', 'III'): 'p10_T50CN3_G2'},
                {**{('C1', 'T50', 'II', m): f'p{11+i}_T50_C1_{m}' for i, m in enumerate(MAINT)},
                 **{('C+', 'T50', 'II', m): f'p{15+i}_T50_Cp_{m}' for i, m in enumerate(MAINT)}}),
   line='SER-02_centreline.gpkg', off=1132.8, srange=(1072.8, 7765.0), half=10.75, frmp='S123',
   crossings=[('EO12 bridge', 434.8), ('parallel road', 1281.4), ('ring road', 1475.4), ('rail', 1544.7),
              ('Serres-Mitrousi road', 1730.3), ('rail Th-A', 2478.8), ('circuit', 3005)]),
 'SER-03': dict(dir='SER-03_AgiosIoannis', name='Agios Ioannis',
   plans=_plans({('T50', 'II'): 'p06_T50_G0', ('T100', 'II'): 'p08_T100_G0', ('T50', 'III'): 'p09_T50CN3_G0'},
                {('T50', 'II'): 'p01_T50_G1', ('T100', 'II'): 'p02_T100_G1', ('T1000', 'II'): 'p03_T1000_G1', ('T50', 'III'): 'p04_T50CN3_G1'},
                {('T50', 'II'): 'p07_T50_G2', ('T100', 'II'): 'p10_T100_G2', ('T50', 'III'): 'p11_T50CN3_G2'},
                {**{('C1', 'T50', 'II', m): f'p{12+i}_T50_C1_{m}' for i, m in enumerate(MAINT)},
                 **{('C+', 'T50', 'II', m): f'p{16+i}_T50_Cp_{m}' for i, m in enumerate(MAINT)},
                 ('C1', 'T50', 'II', 'A13.18km2'): 'p05_T50A13_G1'}),
   line='SER-03_centreline.gpkg', off=0.0, srange=(150, 3252), half=6.5, frmp='S123',
   crossings=[('bridge s236', 226), ('road bridge s271', 261), ('EO12 culvert', 1414), ('bridge s2455', 2445)]),
 'SER-04': dict(dir='SER-04_Krousovitis', name='Krousovitis',
   plans=_plans({('T50', 'II'): 'p06_T50_G0', ('T100', 'II'): 'p08_T100_G0', ('T50', 'III'): 'p09_T50CN3_G0'},
                {('T50', 'II'): 'p01_T50_G1', ('T100', 'II'): 'p02_T100_G1', ('T1000', 'II'): 'p03_T1000_G1', ('T50', 'III'): 'p04_T50CN3_G1'},
                {('T50', 'II'): 'p07_T50_G2', ('T100', 'II'): 'p10_T100_G2', ('T50', 'III'): 'p11_T50CN3_G2'},
                {**{('C1', 'T50', 'II', m): f'p{14+i}_T50_C1_{m}' for i, m in enumerate(MAINT)},
                 **{('C+', 'T50', 'II', m): f'p{18+i}_T50_Cp_{m}' for i, m in enumerate(MAINT)}}),
   line='SER-04_centreline.gpkg', off=2806.5058, srange=None, half=28.0, frmp='S4',
   crossings=[('old EO', 942), ('A25', 1606), ('railway', 1704), ('farm track', 4429), ('Belitsa junction', 10509)]),
}
REF = ('C1', 'T50', 'II', 'base')

def skey(k):
    """readable scenario label"""
    ch, T, cn, m = k
    s = f'{ch} {T}' + (' CN III' if cn == 'III' else '')
    return s + ('' if m == 'base' else f' {m}')

class Case:
    """Grid, domain, centreline stations, reaches and buildings of one case (on the grid of its result rasters)."""
    def __init__(self, c):
        self.c = c; cf = CASES[c]; self.cf = cf
        ref = glob.glob(os.path.join(ROOT, cf['dir'], '03_hecras', 'results', cf['plans'][REF], '*maxdepth.tif'))[0]
        with rasterio.open(ref) as r:
            self.tr, self.sh, self.crs = r.transform, r.shape, r.crs
        self.px = abs(self.tr.a); self.pxa = self.px * self.px / 1e4            # pixel size (m), pixel area (ha)
        per = json.load(open(os.path.join(P1, 'data', 'model_perimeters.json')))[c]['pts']
        self.poly = Polygon(per)
        self.dom = rasterize([(mapping(self.poly), 1)], out_shape=self.sh, transform=self.tr, fill=0, dtype='uint8').astype(bool)
        # centreline stations every 2 m, chainage ch = s - off (case convention)
        L = gpd.read_file(os.path.join(P1, 'data', 'centrelines', cf['line'])).geometry.iloc[0]
        s0, s1 = cf['srange'] if cf['srange'] else (cf['off'] - 1300, L.length - 5)
        ss = np.arange(s0, s1, 2.0); pts = np.array([L.interpolate(s).coords[0] for s in ss])
        self.kt = cKDTree(pts); self.ss = ss - cf['off']
        rr, cc = np.mgrid[0:self.sh[0], 0:self.sh[1]]
        self.X = self.tr.c + self.tr.a * (cc + 0.5); self.Y = self.tr.f + self.tr.e * (rr + 0.5)
        dist, ix = self.kt.query(np.c_[self.X.ravel(), self.Y.ravel()])
        self.dist = dist.reshape(self.sh); self.ch = self.ss[ix].reshape(self.sh)
        # reaches: between crossings, split so that none is longer than 500 m (same rule as D3)
        cut = sorted(x for _, x in cf['crossings']); lo, hi = self.ss[0], self.ss[-1]; e = [lo]
        for x in cut + [hi]:
            while x - e[-1] > 500: e.append(e[-1] + 500)
            if x > e[-1]: e.append(x)
        self.edges = np.array(e)
        self.reach = np.clip(np.digitize(self.ch, self.edges) - 1, 0, len(e) - 2)
        self._bld = None

    def reach_table(self):
        return [dict(case=self.c, reach=f'{self.c}-R{j+1:02d}', ch_from=round(self.edges[j]), ch_to=round(self.edges[j + 1]),
                     length_m=round(self.edges[j + 1] - self.edges[j])) for j in range(len(self.edges) - 1)]

    def reach_lines(self):
        """reach segments of the centreline as a GeoDataFrame (for maps)"""
        from shapely.ops import substring
        L = gpd.read_file(os.path.join(P1, 'data', 'centrelines', self.cf['line'])).geometry.iloc[0]
        rt = self.reach_table(); off = self.cf['off']
        geoms = [substring(L, r['ch_from'] + off, r['ch_to'] + off) for r in rt]
        return gpd.GeoDataFrame(rt, geometry=geoms, crs=self.crs)

    def depth(self, key):
        """max depth (m) of a scenario on the case grid, 0 outside the domain / dry"""
        f = glob.glob(os.path.join(ROOT, self.cf['dir'], '03_hecras', 'results', self.cf['plans'][key], '*maxdepth.tif'))[0]
        with rasterio.open(f) as r:
            if r.transform == self.tr and r.shape == self.sh:
                d = r.read(1).astype('float32'); nd = r.nodata
            else:
                d = np.zeros(self.sh, 'float32'); nd = r.nodata
                reproject(r.read(1).astype('float32'), d, src_transform=r.transform, src_crs=r.crs, dst_transform=self.tr,
                          dst_crs=self.crs, src_nodata=nd, dst_nodata=0, resampling=Resampling.nearest)
        if nd is not None: d[d == nd] = 0
        d[~np.isfinite(d) | (d < 0)] = 0
        d[~self.dom] = 0
        return d

    def frmp(self, T):
        """FRMP EL11 extent (bool) for T in ('0050','0100','1000') on the case grid"""
        f = os.path.join(ROOT, '_data_frmp', 'frmp_el11_extents', f"frmp_{T}_{self.cf['frmp']}.tif")
        a = np.zeros(self.sh, 'uint8')
        with rasterio.open(f) as r:
            reproject(r.read(1).astype('uint8'), a, src_transform=r.transform, src_crs=r.crs, dst_transform=self.tr,
                      dst_crs=self.crs, resampling=Resampling.nearest)
        return (a > 0) & self.dom

    def buildings(self):
        """buildings inside the domain: GeoDataFrame with reach of the nearest centreline station, and a label raster"""
        if self._bld is None:
            b = gpd.read_file(os.path.join(ROOT, '_data_buildings', 'buildings_merged_serres.gpkg'), bbox=self.poly.bounds)
            b = b.to_crs(self.crs) if b.crs != self.crs else b
            b = b[b.representative_point().within(self.poly)].reset_index(drop=True)
            p = b.representative_point(); dist, ix = self.kt.query(np.c_[p.x, p.y])
            b['ch'] = self.ss[ix]; b['dist_m'] = dist
            b['reach_i'] = np.clip(np.digitize(b['ch'], self.edges) - 1, 0, len(self.edges) - 2)
            b['reach'] = [f'{self.c}-R{j+1:02d}' for j in b['reach_i']]
            lab = rasterize([(g, i + 1) for i, g in enumerate(b.geometry)], out_shape=self.sh, transform=self.tr,
                            fill=0, dtype='int32', all_touched=True)
            self._bld = (b, lab)
        return self._bld

    def bld_depth(self, D):
        """max depth over each building footprint (all touched cells)"""
        b, lab = self.buildings(); mx = np.zeros(len(b) + 1, 'float32')
        np.maximum.at(mx, lab.ravel(), D.ravel()); return mx[1:]

    def write(self, path, a, dtype='float32', nodata=None):
        with rasterio.open(path, 'w', driver='GTiff', height=self.sh[0], width=self.sh[1], count=1, dtype=dtype, crs=self.crs,
                           transform=self.tr, nodata=nodata, compress='deflate') as o:
            o.write(a.astype(dtype), 1)
