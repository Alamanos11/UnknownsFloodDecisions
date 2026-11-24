"""D3: stepped-inflow runs -> flooded area, buildings and first out-of-bank step per reach.
usage: d3.py CASE"""
import numpy as np, rasterio, json, glob, sys, os, csv, geopandas as gpd
from rasterio.features import rasterize
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, mapping
import scen_post as SP   # reuses grid/terrain/bld_max; CFG
U=SP.U; SL=SP.SL
F=np.array([0.05,0.075,0.10,0.125,0.15,0.175,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.75,1.00])
Q={'SER-01':(149.16,204.13,470.79),'SER-02':(67.0,90.99,207.81),'SER-03':(19.96,26.26,55.95),'SER-04':(428.99,579.24,1306.45)}  # HMS outlet T50, T100, T1000 (CN II)
D=dict(
 **{'SER-01':dict(line=os.path.join(SP.U, 'SER-P1_paper_decisions', 'data', 'centrelines', 'SER-01_centreline.gpkg'),off=0.0,sr=(40,5084),X=[('K1 Meg. Alexandrou bridge',2055),('K2 road',3030),('K3 road',3255),('K4 railway',3875)],half=9.5,plans={'C0':'p21','C1':'p22','Cp':'p23'},c0=['Terrain/DEM_open.T0_open_1m.tif']),
    'SER-02':dict(line=os.path.join(SP.U, 'SER-P1_paper_decisions', 'data', 'centrelines', 'SER-02_centreline.gpkg'),off=1132.8,sr=(1072.8,7765.0),X=[('EO12 bridge',434.8),('parallel road',1281.4),('ring road',1475.4),('rail Zachareos',1544.7),('Serres-Mitrousi road',1730.3),('rail Th-A',2478.8),('circuit',3005)],half=10.75,plans={'C0':'p19','C1':'p20','Cp':'p21'},c0='Terrain/DEM_open.T0_DEMopen_2m.tif'),
    'SER-03':dict(line=os.path.join(SP.U, 'SER-P1_paper_decisions', 'data', 'centrelines', 'SER-03_centreline.gpkg'),off=0.0,sr=(150,3252),X=[('bridge s236',226),('road bridge s271',261),('EO12 culvert',1414),('bridge s2455',2445)],half=6.5,plans={'C0':'p20','C1':'p21','Cp':'p22'},c0='Terrain/DEM_open.T0_DEMopen_2m.tif'),
    'SER-04':dict(line=os.path.join(SP.U, 'SER-P1_paper_decisions', 'data', 'centrelines', 'SER-04_centreline.gpkg'),off=2806.5058,sr=None,X=[('old EO',942),('A25',1606),('railway',1704),('farm track',4429),('Belitsa junction',10509)],half=28,plans={'C0':'p22','C1':'p23','Cp':'p24'},c0='geometry/T0_DEM_5m.tif')})
c=sys.argv[1]; d=D[c]; cf=SP.CFG[c]; tr,sh,crs,ref=SP.grid(c); px=abs(tr.a*tr.e)/1e4
poly=Polygon(SP.P[c]['pts']); dom=rasterize([(mapping(poly),1)],out_shape=sh,transform=tr,fill=0,dtype='uint8').astype(bool)
L=gpd.read_file(d['line']).geometry[0]; s0,s1=d['sr'] if d['sr'] else (d['off']-1300,L.length-5)
ss=np.arange(s0,s1,2.0); pts=np.array([L.interpolate(s).coords[0] for s in ss]); kt=cKDTree(pts)
rows_,cols_=np.mgrid[0:sh[0],0:sh[1]]; XX=tr.c+tr.a*(cols_+0.5); YY=tr.f+tr.e*(rows_+0.5)
dist,ix=kt.query(np.c_[XX.ravel(),YY.ravel()]); dist=dist.reshape(sh); sta=(ss[ix]-d['off']).reshape(sh)   # chainage in the case's own convention
corr=dist<=d['half']; near=(dist<=150)&~corr&dom
# reaches: between crossings, split to <= 500 m
cut=sorted([x for _,x in d['X']]); lo,hi=s0-d['off'],s1-d['off']; edges=[lo]
for x in cut+[hi]:
    while x-edges[-1]>500: edges.append(edges[-1]+500)
    if x>edges[-1]: edges.append(x)
edges=np.array(edges); reach=np.digitize(sta,edges)-1
geo=f"{U}/{cf['dir']}/03_hecras"
def zC0():
    fs=d['c0'] if isinstance(d['c0'],list) else [d['c0']]
    return SP.terrain(c,f"{geo}/{fs[0]}",tr,sh,crs)
Z={'C1':SP.terrain(c,f"{geo}/geometry/{cf['tif']['C1']}.tif",tr,sh,crs),'Cp':SP.terrain(c,f"{geo}/geometry/{cf['tif']['Cp']}.tif",tr,sh,crs)}
rows=[]; first=[]
os.makedirs(f'{c}/D3',exist_ok=True)
for ch,p in d['plans'].items():
    f=f"{SL}/{c}_{p}.npz"
    if not os.path.exists(f): print('missing',f); continue
    if ch=='C0':
        z=zC0()
        if z is None: print('C0 terrain not set for',c); continue
    else: z=Z[ch]
    n=np.load(f); a=[k[:-4] for k in n.files if k.endswith('__cc')][0]; cc=n[a+'__cc']; SW=n[a+'__step_wmax']
    t=cKDTree(cc); dd,ii=t.query(np.c_[XX.ravel(),YY.ravel()],distance_upper_bound=cf['ub']); ok=np.isfinite(dd).reshape(sh); ii=np.where(np.isfinite(dd),ii,0).reshape(sh)
    firstk={}; stack=np.zeros(sh,'int8')+99
    for k in range(15):
        Dk=np.where(ok,SW[k][ii]-z,0); Dk[~(Dk>0)]=0; Dk=Dk.astype('float32')
        w=(Dk>=0.1)&dom; B=SP.bld_max(c,Dk,tr,sh,dom)
        stack[(stack==99)&w]=k
        oob=w&near
        r=dict(case=c,channel=ch,plan=p,step=k+1,frac_T1000=F[k],Q_outlet=round(F[k]*Q[c][2],1),A_wet_ha=round(w.sum()*px,1),A_oob150_ha=round(oob.sum()*px,2),
               bld010=int((B>=0.1).sum()),bld050=int((B>=0.5).sum()))
        rows.append(r)
        for j in range(len(edges)-1):
            if j in firstk: continue
            m=oob&(reach==j)
            if m.sum()*px>=0.25: firstk[j]=k
    with rasterio.open(f'{c}/D3/{p}_Step_{ch}_first_wet_step.tif','w',driver='GTiff',height=sh[0],width=sh[1],count=1,dtype='int8',crs=crs,transform=tr,nodata=99,compress='deflate') as o: o.write(np.where(dom,stack,99).astype('int8'),1)
    for j in range(len(edges)-1):
        k=firstk.get(j)
        first.append(dict(case=c,channel=ch,reach=j+1,ch_from=round(edges[j]),ch_to=round(edges[j+1]),first_oob_step=(k+1 if k is not None else ''),
             frac_T1000=(F[k] if k is not None else ''),Q_outlet=(round(F[k]*Q[c][2],1) if k is not None else ''),
             Q_lower=(round(F[k-1]*Q[c][2],1) if k else '')))
for nm,R in (('D3_steps',rows),('D3_first_overflow',first)):
    if R:
        with open(f'{c}/D3/{c}_{nm}.csv','w',newline='') as fo: w_=csv.DictWriter(fo,fieldnames=list(R[0].keys())); w_.writeheader(); w_.writerows(R)
print('T50/T100 as fraction of T1000:',[round(x/Q[c][2],3) for x in Q[c][:2]])
for r in rows: print(r['channel'],r['step'],r['frac_T1000'],r['Q_outlet'],'A',r['A_wet_ha'],'oob',r['A_oob150_ha'],'bld',r['bld010'],r['bld050'])
for r in first: print(r)
