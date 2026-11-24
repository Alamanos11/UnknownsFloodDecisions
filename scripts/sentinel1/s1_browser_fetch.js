// Sentinel-1 extraction for SER-04, run in the browser console of the Copernicus Browser
// (https://browser.dataspace.copernicus.eu), logged in, with the Sentinel-1 IW VV+VH collection and
// the "VV - decibel gamma0" layer selected (orthorectified, gamma0 ellipsoid, Copernicus DEM, no speckle filter).
// It reuses the page's own Sentinel Hub layer object (S1GRDCDASLayer) and its own request config,
// so no credentials are handled by this script. Run on 23/11/2026 (Copernicus Browser v1.65.0).
//
// Steps used:
// 1. Find the Leaflet map through the React fiber of .leaflet-container and the Sentinel Hub tile layer
//    (options.tileSize 512, has .layer promise) -> window.__L, window.__U = await window.__L.layer.
// 2. Wrap __U.getMap once to capture one real tile request: args [params, apiType, reqConfig] -> window.__cap[0].
// 3. Request the SER-04 box (EPSG:3857, lon 23.285-23.410, lat 41.155-41.260) at 13.3 m Mercator
//    (~10 m ground) with the evalscript below, one request per day and orbit direction.

const R = 6378137, mx = lon => R * lon * Math.PI / 180, my = lat => R * Math.log(Math.tan(Math.PI / 4 + lat * Math.PI / 360));
const x0 = mx(23.285), x1 = mx(23.410), y0 = my(41.155), y1 = my(41.260), px = 13.3;
const W = Math.round((x1 - x0) / px), H = Math.round((y1 - y0) / px);   // 1046 x 1168

const EV = `//VERSION=3
function setup(){return {input:["VV","dataMask"],output:{bands:2,sampleType:"UINT8"}};}
function evaluatePixel(s){ if(!s.dataMask||s.VV<=0) return [0,0]; let db=10*Math.log10(s.VV);
  let c=Math.round((db+35)*6); c=Math.max(1,Math.min(255,c)); return [c,255]; }`;   // code = (dB+35)*6

async function fetchDB(day, orbit) {          // orbit 'DESCENDING' | 'ASCENDING'
  const u = window.__U, base = window.__cap[0], BB = base[0].bbox.constructor;
  const p = { ...base[0], bbox: new BB(base[0].bbox.crs, x0, y0, x1, y1), width: W, height: H, format: 'image/png',
              fromTime: new Date(day + 'T00:00:00Z'), toTime: new Date(day + 'T23:59:59Z'), orbitDirection: orbit,
              preview: undefined, effects: undefined };
  const ev = u.evalscript, od = u.orbitDirection; u.evalscript = EV; u.orbitDirection = orbit;
  let blob; try { blob = await u.__orig.call(u, p, base[1], base[2]); } finally { u.evalscript = ev; u.orbitDirection = od; }
  const bm = await createImageBitmap(blob), cv = new OffscreenCanvas(bm.width, bm.height), cx = cv.getContext('2d');
  cx.drawImage(bm, 0, 0); const d = cx.getImageData(0, 0, bm.width, bm.height).data, a = new Uint8Array(W * H);
  for (let i = 0, j = 0; j < a.length; i += 4, j++) a[j] = d[i + 3] > 0 ? d[i] : 0;
  return a;
}

// 4. Per image: dB = code/6 - 35; 3x3 boxcar in linear power.
// 5. References (pre-event, same track): D21 = mean(27/12/2020, 2/1/2021 desc.), A21 = 2/1/2021 asc.,
//    D26 = mean(12/1, 18/1/2026 desc.), A26 = mean(12/1, 18/1/2026 asc.).
// 6. New open water = post <= -18 dB AND (ref - post) >= 3 dB.  refdark_D21 = ref D21 <= -18 dB (permanent dark).
// 7. Masks bit-packed, gzip, base64 -> s1_masks.json (meta holds grid W,H,x0..y1, EPSG:3857).
// The 16:15 UTC ascending track (28/12, 3/1, 9/1, 15/1, 21/1) covers only ~2 % of the box and was not used.
