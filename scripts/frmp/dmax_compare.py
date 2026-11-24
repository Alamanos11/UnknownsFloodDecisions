import numpy as np, rasterio, json, glob, re, csv
from rasterio.warp import reproject, Resampling
from rasterio.features import rasterize
from shapely.geometry import Polygon, mapping
import os; U=os.environ.get("P1_ROOT", "02_CASES")   # set P1_ROOT to the project folder 02_CASES
P=json.load(open(f"{U}/_data_frmp/model_perimeters.json"))
CASES={'SER-01':('SER-01_AgioiAnargyroi','S123',{'G2':'C1','G1':'C+','G0open':'C0'}),
       'SER-02':('SER-02_AgiaVarvara','S123',{'G1':'C1','G2':'C+','G0':'C0'}),
       'SER-03':('SER-03_AgiosIoannis','S123',{'G1':'C1','G2':'C+','G0':'C0'}),
       'SER-04':('SER-04_Krousovitis','S4',{'G1':'C1','G2':'C+','G0':'C0'})}
TT={'T50':'0050','T100':'0100','T1000':'1000'}
edges=[0.1,0.2,0.5,1.0,2.0,1e9]   # model classes start at 0.10 (our wet threshold)
def cls(d,lo):
    c=np.zeros(d.shape,np.int8)  # 0 dry
    for i,(a,b) in enumerate(zip([lo]+edges[1:-1],edges[1:])): c[(d>=a)&(d<b)]=i+1
    return c
rows=[]
for c,(folder,box,gm) in CASES.items():
    poly=Polygon(P[c]['pts'])
    for f in sorted(glob.glob(f"{U}/{folder}/03_hecras/results/p*/*maxdepth.tif")):
        plan=f.split('/')[-2]; m=re.match(r'p\d+_(\w+?)_(G\w+)$',plan)
        if not m or m.group(1) not in TT: continue
        run,ch=m.group(1),gm.get(m.group(2)); T=TT[run]
        with rasterio.open(f) as d: dep=d.read(1).astype('float32'); tr=d.transform; sh=dep.shape; crs=d.crs; nd=d.nodata
        if nd is not None: dep[dep==nd]=0
        dep[~np.isfinite(dep)]=0
        dom=rasterize([(mapping(poly),1)],out_shape=sh,transform=tr,fill=0,dtype='uint8').astype(bool)
        fr=np.full(sh,-1,'float32')
        with rasterio.open(f'dmax/frmp_dmax_{T}_{box}.tif') as s:
            reproject(s.read(1),fr,src_transform=s.transform,src_crs=s.crs,dst_transform=tr,dst_crs=crs,resampling=Resampling.nearest,src_nodata=-1,dst_nodata=-1)
        px=abs(tr.a*tr.e)/1e4
        mw=(dep>=0.10)&dom; fw=(fr>=0)&dom
        both=mw&fw
        cm=cls(np.where(mw,dep,0),0.10); cf=np.where(fw,cls(np.maximum(fr,0.1),0.10),0)  # FRMP <0.2 class merged as 'shallow'
        # 4 hazard classes: <0.5, 0.5-1, 1-2, >=2 for both wet
        def h(x): return np.digitize(x,[0.5,1.0,2.0])
        hm=h(dep[both]); hf=h(fr[both])
        dd=np.minimum(dep[both],2.0)-fr[both]
        r=dict(case=c,plan=plan,run=run,channel=ch,frmp_T=int(T),
               A_model_ha=round(mw.sum()*px,1),A_frmp_dmax_ha=round(fw.sum()*px,1),A_both_ha=round(both.sum()*px,1),
               model_ge05_ha=round(((dep>=0.5)&dom).sum()*px,1),frmp_ge05_ha=round(((fr>=0.5)&dom).sum()*px,1),
               model_ge1_ha=round(((dep>=1)&dom).sum()*px,1),frmp_ge1_ha=round(((fr>=1)&dom).sum()*px,1),
               same_class_share=round(float((hm==hf).mean()),3) if both.any() else None,
               within1class_share=round(float((abs(hm-hf)<=1).mean()),3) if both.any() else None,
               dd_median_m=round(float(np.median(dd)),2) if both.any() else None,
               dd_p10_m=round(float(np.percentile(dd,10)),2) if both.any() else None,
               dd_p90_m=round(float(np.percentile(dd,90)),2) if both.any() else None,
               model_deeper_share=round(float((dd>0.1).mean()),3) if both.any() else None,
               frmp_deeper_share=round(float((dd<-0.1).mean()),3) if both.any() else None)
        rows.append(r); print(r)
with open('out/frmp_dmax_vs_model.csv','w',newline='') as fo:
    w=csv.DictWriter(fo,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
