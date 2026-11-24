// Run in a browser tab on https://gis.floods.ypeka.gr (same origin). Fetches the EL11 max-depth layer
// with the server's own style EL_DMAX (colour ramp 0 / 0.4 / 1 / 1.5 / 2 m), inverts the ramp per pixel,
// and stores gzip+base64 of a uint8 grid: 0 = no data (dry), 1..201 = depth (cm)+1 (0-2.00 m), 250 = >= 2 m.
// SLD_BODY requests (custom styles) are refused by the server's firewall (403), and WCS is off,
// so the published style is the only way to read values. Check: 0 pixels off the ramp in all six grids.
async function getMap(ly,bbox,w,h){const u='https://gis.floods.ypeka.gr/geoserver/floods_rev1/wms?service=WMS&version=1.1.1&request=GetMap&layers='+encodeURIComponent(ly)+'&styles=EL_DMAX&srs=EPSG:2100&bbox='+bbox+'&width='+w+'&height='+h+'&format=image/png&transparent=true&interpolations=nearest%20neighbor&format_options=antialiasing:off'; const r=await fetch(u); return await r.blob();}
async function decode(ly,box){ // ly e.g. 'EL11_DMΑΧ_0100' (Greek Α and Χ in the layer name); box [x0,y0,x1,y1] EGSA87, 5 m pixels
 const [x0,y0,x1,y1]=box, W=Math.round((x1-x0)/5), H=Math.round((y1-y0)/5);
 const im=await createImageBitmap(await getMap(ly,box.join(','),W,H)); const c=new OffscreenCanvas(W,H); const g=c.getContext('2d'); g.drawImage(im,0,0);
 const d=g.getImageData(0,0,W,H).data, q=[0,0.4,1,1.5,2], C=[[182,237,240],[116,180,232],[31,131,224],[29,68,184],[9,9,145]], o=new Uint8Array(W*H);
 for(let i=0,j=0;i<d.length;i+=4,j++){ if(d[i+3]==0) continue; const p=[d[i],d[i+1],d[i+2]]; let best=1e9,dep=0;
  for(let k=0;k<4;k++){const a=C[k],b=C[k+1],D=[b[0]-a[0],b[1]-a[1],b[2]-a[2]],L=D[0]**2+D[1]**2+D[2]**2; let t=((p[0]-a[0])*D[0]+(p[1]-a[1])*D[1]+(p[2]-a[2])*D[2])/L; t=Math.max(0,Math.min(1,t));
   const e=[p[0]-a[0]-t*D[0],p[1]-a[1]-t*D[1],p[2]-a[2]-t*D[2]], dd=e[0]**2+e[1]**2+e[2]**2; if(dd<best){best=dd; dep=q[k]+t*(q[k+1]-q[k]);}}
  o[j]=(p[0]==9&&p[1]==9&&p[2]==145)?250:1+Math.min(200,Math.round(dep*100)); }
 const ab=new Uint8Array(await new Response(new Blob([o]).stream().pipeThrough(new CompressionStream('gzip'))).arrayBuffer());
 let s=''; for(let i=0;i<ab.length;i+=0x8000) s+=String.fromCharCode.apply(null,ab.subarray(i,i+0x8000)); return btoa(s); }
