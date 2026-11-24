import numpy as np, rasterio, json, glob, os, re, csv
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from rasterio.warp import reproject, Resampling
from rasterio.features import rasterize
from shapely.geometry import Polygon, mapping
import os; U=os.environ.get("P1_ROOT", "02_CASES")   # set P1_ROOT to the project folder 02_CASES
P=json.load(open(f"{U}/_data_frmp/model_perimeters.json"))
CASES={'SER-01':('SER-01_AgioiAnargyroi','S123',{'G2':'C1','G1':'C+','G0open':'C0'}),
       'SER-02':('SER-02_AgiaVarvara','S123',{'G1':'C1','G2':'C+','G0':'C0'}),
       'SER-03':('SER-03_AgiosIoannis','S123',{'G1':'C1','G2':'C+','G0':'C0'}),
       'SER-04':('SER-04_Krousovitis','S4',{'G1':'C1','G2':'C+','G0':'C0'})}
rows=[]
def frmp_on(tr,sh,crs,T,box):
    fr=np.zeros(sh,'uint8')
    with rasterio.open(f'frmp_{T}_{box}.tif') as s:
        reproject(s.read(1),fr,src_transform=s.transform,src_crs=s.crs,dst_transform=tr,dst_crs=crs,resampling=Resampling.nearest)
    return fr.astype(bool)
for c,(folder,box,gm) in CASES.items():
    poly=Polygon(P[c]['pts'])
    for f in sorted(glob.glob(f"{U}/{folder}/03_hecras/results/p*/*maxdepth.tif")):
        plan=f.split('/')[-2]; m=re.match(r'p\d+_(\w+?)_(G\w+)$',plan)
        if not m: continue
        run,g=m.group(1),m.group(2); ch=gm.get(g,g)
        T={'T50':'0050','T100':'0100','T1000':'1000','T50CN3':'0050','T50A13':'0050'}.get(run)
        if T is None: continue
        with rasterio.open(f) as d: dep=d.read(1); tr=d.transform; sh=dep.shape; crs=d.crs
        dom=rasterize([(mapping(poly),1)],out_shape=sh,transform=tr,fill=0,dtype='uint8').astype(bool)
        ours=(dep>=0.10)&dom; fr=frmp_on(tr,sh,crs,T,box)&dom
        px=abs(tr.a*tr.e)/1e4
        A=ours.sum()*px; B=fr.sum()*px; I=(ours&fr).sum()*px; Un=(ours|fr).sum()*px
        rows.append(dict(case=c,plan=plan,run=run,channel=ch,frmp_T=int(T),A_model_ha=round(A,1),A_frmp_ha=round(B,1),overlap_ha=round(I,1),
                         CSI=round(I/Un,3) if Un else None, hit_rate=round(I/B,3) if B else None, model_in_frmp=round(I/A,3) if A else None, bias=round(A/B,2) if B else None))
with open('out/frmp_vs_model.csv','w',newline='') as fo:
    w=csv.DictWriter(fo,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
