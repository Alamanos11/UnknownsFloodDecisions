"""Map maintenance-scenario runs (slim npz) to max-depth rasters and compare with the base T50 run of the same channel.
usage: scen_post.py [CASE ...]"""
import numpy as np, rasterio, json, glob, sys, os, csv, geopandas as gpd
from rasterio.warp import reproject, Resampling
from rasterio.features import rasterize
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, mapping
import os; U=os.environ.get("P1_ROOT", "02_CASES")   # set P1_ROOT to the project folder 02_CASES; SL=f"{U}/_scen_build/slim"; TO=os.environ.get("SCEN_TERRAIN", f"{U}/_scen_build/terrain")   # 1-5 m channel terrains of C1/C+ (HEC-RAS RAS Mapper exports)
P=json.load(open(f"{U}/_data_frmp/model_perimeters.json"))
CFG={'SER-01':dict(dir='SER-01_AgioiAnargyroi',tif={'C1':'T2_channelB_1m','Cp':'T1_channel_1m'},base5='T0_dem5m_clip',ub=12,base={'C1':'p01','Cp':'p10'},first=13),
     'SER-02':dict(dir='SER-02_AgiaVarvara',tif={'C1':'T1_DEMchannel_2m_v2','Cp':'T2_DEMchannelB_2m'},base5=None,ub=12,base={'C1':'p01','Cp':'p08'},first=11),
     'SER-03':dict(dir='SER-03_AgiosIoannis',tif={'C1':'T1_DEMchannel_2m_v2','Cp':'T2_DEMchannelB_2m_v2'},base5=None,ub=12,base={'C1':'p01','Cp':'p07'},first=12),
     'SER-04':dict(dir='SER-04_Krousovitis',tif={'C1':'T1_DEMchannel_5m','Cp':'T2_DEMchannelB_5m'},base5=None,ub=45,base={'C1':'p01','Cp':'p07'},first=14)}
SC=[('C1','VegOver'),('C1','VegClear'),('C1','Block50'),('C1','Sed05'),('Cp','VegOver'),('Cp','VegClear'),('Cp','Block50'),('Cp','Sed05')]
BLD={c:gpd.read_file(os.path.join(os.environ.get("BLD_DIR", f"{U}/_data_buildings"), f"bld_merged_{c}.gpkg")) for c in CFG}   # buildings clipped per case
def grid(c):
    f=sorted(glob.glob(f"{U}/{CFG[c]['dir']}/03_hecras/results/p01*/*maxdepth.tif"))[0]
    r=rasterio.open(f); return r.transform,r.shape,r.crs,f
def terrain(c,src,tr,sh,crs):
    z=np.full(sh,np.nan,'float32')
    def put(fn,rs):
        with rasterio.open(fn) as s:
            a=np.full(sh,np.nan,'float32'); reproject(s.read(1).astype('float32'),a,src_transform=s.transform,src_crs=s.crs,src_nodata=s.nodata,dst_transform=tr,dst_crs=crs,dst_nodata=np.nan,resampling=rs)
        return a
    if CFG[c]['base5']: z=put(f"{U}/{CFG[c]['dir']}/03_hecras/geometry/{CFG[c]['base5']}.tif",Resampling.bilinear)
    a=put(src,Resampling.average); z=np.where(np.isfinite(a),a,z); return z
def mapdepth(c,npz,z,tr,sh):
    d=np.load(npz); a=[k[:-4] for k in d.files if k.endswith('__cc')][0]
    cc=d[a+'__cc']; w=d[a+'__wmax']; t=cKDTree(cc)
    rows,cols=np.mgrid[0:sh[0],0:sh[1]]; X=tr.c+tr.a*(cols+0.5); Y=tr.f+tr.e*(rows+0.5)
    dd,ii=t.query(np.c_[X.ravel(),Y.ravel()],distance_upper_bound=CFG[c]['ub'])
    ok=np.isfinite(dd).reshape(sh); ii=np.where(np.isfinite(dd),ii,0).reshape(sh)
    D=np.where(ok,w[ii]-z,np.nan); D[~(D>0)]=0; return D.astype('float32'), json.loads(str(d['meta']))
def bld_max(c,D,tr,sh,dom):
    b=BLD[c]; lab=rasterize([(g,i+1) for i,g in enumerate(b.geometry)],out_shape=sh,transform=tr,fill=0,dtype='int32',all_touched=True)
    mx=np.zeros(len(b)+1); np.maximum.at(mx,lab.ravel(),np.where(dom,D,0).ravel()); return mx[1:]
if __name__=='__main__':
    rows=[]
    for c in (sys.argv[1:] or list(CFG)):
        cf=CFG[c]; tr,sh,crs,ref=grid(c); px=abs(tr.a*tr.e)/1e4; os.makedirs(c,exist_ok=True)
        poly=Polygon(P[c]['pts']); dom=rasterize([(mapping(poly),1)],out_shape=sh,transform=tr,fill=0,dtype='uint8').astype(bool)
        edge=dom & ~rasterize([(mapping(poly.buffer(-30)),1)],out_shape=sh,transform=tr,fill=0,dtype='uint8').astype(bool)
        geo=f"{U}/{cf['dir']}/03_hecras/geometry"
        Z={ch:terrain(c,f"{geo}/{cf['tif'][ch]}.tif",tr,sh,crs) for ch in ('C1','Cp')}
        base={}
        for ch in ('C1','Cp'):
            D,m=mapdepth(c,f"{SL}/{c}_{cf['base'][ch]}.npz",Z[ch],tr,sh); base[ch]=(D,bld_max(c,D,tr,sh,dom))
            if ch=='C1':  # check against the delivered core map
                with rasterio.open(ref) as r: R=r.read(1); R=np.where(R>0,R,0)
                print(c,'mapping check p01: mine %.1f ha, delivered %.1f ha, |diff| cells %.3f'%(((D>=0.1)&dom).sum()*px,((R>=0.1)&dom).sum()*px,(((D>=0.1)!=(R>=0.1))&dom).mean()/max(1e-9,((R>=0.1)&dom).mean())))
        for k,(ch,sc) in enumerate(SC):
            p=f"p{cf['first']+k:02d}"; f=f"{SL}/{c}_{p}.npz"
            if not os.path.exists(f): continue
            zz=Z[ch] if sc in ('VegOver','VegClear') else terrain(c,f"{TO}/{c.replace('-','')}_{ch}_{'BK' if sc=='Block50' else 'SD'}.tif",tr,sh,crs)
            D,m=mapdepth(c,f,zz,tr,sh); Db,Bb=base[ch]; Bm=bld_max(c,D,tr,sh,dom)
            with rasterio.open(f"{c}/{p}_T50_{ch}_{sc}_maxdepth.tif",'w',driver='GTiff',height=sh[0],width=sh[1],count=1,dtype='float32',crs=crs,transform=tr,nodata=-9999,compress='deflate',predictor=3) as o: o.write(np.where(dom,D,-9999).astype('float32'),1)
            w=(D>=0.1)&dom; wb=(Db>=0.1)&dom
            rows.append(dict(case=c,plan=p,channel=ch,scenario=sc,solution=m['summary'].get('Solution'),run_time=m['summary'].get('Computation Time Total'),
                A_ha=round(w.sum()*px,1),A_base_ha=round(wb.sum()*px,1),dA_ha=round((w.sum()-wb.sum())*px,1),dA_pct=round(100*(w.sum()/wb.sum()-1),1),
                new_wet_ha=round((w&~wb).sum()*px,1),new_dry_ha=round((wb&~w).sum()*px,1),
                A05_ha=round(((D>=0.5)&dom).sum()*px,1),A05_base_ha=round(((Db>=0.5)&dom).sum()*px,1),
                bld010=int((Bm>=0.1).sum()),bld010_base=int((Bb>=0.1).sum()),bld050=int((Bm>=0.5).sum()),bld050_base=int((Bb>=0.5).sum()),
                bld_new010=int(((Bm>=0.1)&(Bb<0.1)).sum()),bld_lost010=int(((Bm<0.1)&(Bb>=0.1)).sum()),
                edge30_wet_ha=round((w&edge).sum()*px,2),edge30_wet_base_ha=round((wb&edge).sum()*px,2)))
            print(rows[-1])
    if rows:
        fn='scen_D2_summary.csv'; old=[r for r in csv.DictReader(open(fn))] if os.path.exists(fn) else []
        keep=[r for r in old if (r['case'],r['plan']) not in {(x['case'],x['plan']) for x in rows}]
        allr=keep+[{k:str(v) for k,v in r.items()} for r in rows]; allr.sort(key=lambda r:(r['case'],r['plan']))
        with open(fn,'w',newline='') as fo: w_=csv.DictWriter(fo,fieldnames=list(rows[0].keys())); w_.writeheader(); w_.writerows(allr)
