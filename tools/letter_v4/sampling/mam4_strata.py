#!/usr/bin/env python3
"""Stratified sampling of MAM4 surface states: 14 regions x 4 seasons x {typical, polluted (P90)}.

Usage: python mam4_strata.py <hist_dir> <case> <regimes.npz> <landscan_2010.tif> <out_dir>

typical  = the (cell, day) of the region-season closest to the region-season median feature vector
polluted = among (cell, day) with the region's pollution metric >= its region-season P90, the one
           closest to the median feature vector of that subset. Metric: fine PM (a1+a2+a4 dry mass);
           total dust mass (incl. coarse) for the dust regions.
Weights per stratum: area x time, LandScan-2010 population x time, low-cloud x area x time.
Also writes population weights of the 30 area-based regimes (from their per-cell counts).
"""
import sys, os, glob, json, time
import numpy as np
from netCDF4 import Dataset
import rasterio

HIST, CASE, REG, POP, OUT = sys.argv[1:6]
os.makedirs(OUT, exist_ok=True)
RG = np.load(REG, allow_pickle=True)
RAW = list(RG['raw_names']); names = list(RG['names']); mu, sd, wfeat = RG['mu'], RG['sd'], RG['wfeat']
lat, lon = RG['lat'], RG['lon']; nlat, nlon = len(lat), len(lon); ncell = nlat * nlon
LON2D, LAT2D = np.meshgrid(lon, lat); LONf, LATf = LON2D.ravel(), LAT2D.ravel()
with Dataset(sorted(glob.glob(f'{HIST}/{CASE}.cam.h0.0002-*.nc'))[0]) as nc:
    gw = np.asarray(nc['gw'][:])
AREA = np.repeat(gw, nlon)

# ---------------- population on the f19 grid (sum of LandScan 1 km counts per cell) ---------------
pop_file = f'{OUT}/pop_f19.npy'
if os.path.exists(pop_file):
    pop = np.load(pop_file)
else:
    pop = np.zeros(ncell)
    latb = np.concatenate([[-90], 0.5 * (lat[1:] + lat[:-1]), [90]])
    dlon = lon[1] - lon[0]
    with rasterio.open(POP) as ds:
        nodata = ds.nodata; T = ds.transform; H, Wd = ds.height, ds.width
        xs = T.c + (np.arange(Wd) + 0.5) * T.a
        ilon = np.round(np.mod(xs, 360) / dlon).astype(int) % nlon
        for r0 in range(0, H, 1200):
            win = rasterio.windows.Window(0, r0, Wd, min(1200, H - r0))
            a = ds.read(1, window=win).astype(np.float64)
            a[(a < 0) | (a == nodata)] = 0
            ys = T.f + (r0 + np.arange(a.shape[0]) + 0.5) * T.e
            ilat = np.clip(np.searchsorted(latb, ys) - 1, 0, nlat - 1)
            for i in np.unique(ilat):
                rows = a[ilat == i].sum(0)
                np.add.at(pop, i * nlon + ilon, rows)
    np.save(pop_file, pop)
print('population total (billion)', pop.sum() / 1e9, flush=True)

# ---------------- regions -------------------------------------------------------------------------
def box(la0, la1, lo0, lo1):
    lo = np.mod(LONf, 360); lo0, lo1 = lo0 % 360, lo1 % 360
    inlon = (lo >= lo0) & (lo <= lo1) if lo0 <= lo1 else (lo >= lo0) | (lo <= lo1)
    return (LATf >= la0) & (LATf <= la1) & inlon
REGIONS = {  # name: (mask, metric, kind)
    'north_china_plain_east_china': (box(28, 42, 110, 122), 'pm', 'land'),
    'indo_gangetic_plain': (box(22, 32, 72, 90), 'pm', 'land'),
    'southeast_asia': (box(8, 25, 95, 110), 'pm', 'land'),
    'western_europe': (box(43, 55, -10, 15), 'pm', 'land'),
    'eastern_us': (box(32, 45, -90, -70), 'pm', 'land'),
    'amazon': (box(-15, 0, -65, -45), 'pm', 'land'),
    'central_africa_biomass': (box(-15, 5, 10, 35), 'pm', 'land'),
    'sahara_sahel_dust': (box(12, 30, -15, 30), 'dust', 'land'),
    'east_asian_dust': (box(36, 45, 75, 110), 'dust', 'land'),
    'middle_east_dust': (box(20, 35, 40, 60), 'dust', 'land'),
    'nw_pacific_outflow': (box(25, 45, 125, 160), 'pm', 'ocean'),
    'north_atlantic': (box(40, 60, -40, -10), 'pm', 'ocean'),
    'southern_ocean': (box(-65, -45, 0, 359.99), 'pm', 'ocean'),
    'arctic': (box(70, 90, 0, 359.99), 'pm', 'any'),
    'southeast_atlantic_stratocumulus': (box(-20, 0, -15, 15), 'pm', 'ocean')}
SEAS = {'DJF': (12, 1, 2), 'MAM': (3, 4, 5), 'JJA': (6, 7, 8), 'SON': (9, 10, 11)}
cells = np.flatnonzero(np.any([m for m, _, _ in REGIONS.values()], axis=0))
cpos = {c: i for i, c in enumerate(cells)}
print('region cells', len(cells), {k: int(m.sum()) for k, (m, _, _) in REGIONS.items()}, flush=True)

# ---------------- read daily surface states for region cells (years 2-5, every day) --------------
files = []
for y in range(2, 6):
    files += sorted(glob.glob(f'{HIST}/{CASE}.cam.h0.{y:04d}-*.nc'))
nd = len(files)
raw = np.empty((nd, len(cells), len(RAW)), np.float32); months = np.empty(nd, int); dates = []
t0 = time.time()
for i, f in enumerate(files):
    with Dataset(f) as nc:
        for j, v in enumerate(RAW):
            x = nc[v]; a = np.asarray(x[0, -1] if x.ndim == 4 else x[0], dtype=np.float32).ravel()
            raw[i, :, j] = a[cells]
    d = os.path.basename(f).split('.h0.')[1][:10]; dates.append(d); months[i] = int(d[5:7])
    if i % 200 == 0:
        print(i, f'{time.time() - t0:.0f}s', flush=True)
ix = {v: j for j, v in enumerate(RAW)}
R = raw.reshape(nd * len(cells), len(RAW))
rho = R[:, ix['PS']] / (287.05 * R[:, ix['TREFHT']])
ug = lambda vs: sum(R[:, ix[v + '_SRF']] for v in vs) * rho * 1e9
SPEC = {'BC': ['bc_a1', 'bc_a4'], 'POM': ['pom_a1', 'pom_a4'], 'SOA': [f'soa{i}_a{m}' for i in range(1, 6) for m in (1, 2)],
        'SO4': ['so4_a1', 'so4_a2', 'so4_a3'], 'DST': ['dst_a1', 'dst_a2', 'dst_a3'], 'NCL': ['ncl_a1', 'ncl_a2', 'ncl_a3']}
NUM = ['num_a1', 'num_a2', 'num_a3', 'num_a4']
fine = {k: ug([v for v in vs if not v.endswith('_a3')]) for k, vs in SPEC.items()}
coarse = sum(ug([v for v in vs if v.endswith('_a3')]) for vs in SPEC.values())
mfine = sum(fine.values()) + 1e-6
dust_tot = ug(SPEC['DST'])
ntot = sum(R[:, ix[n + '_SRF']] for n in NUM) * rho * 1e-6
ppb = lambda g: R[:, ix[g + '_SRF']] * 1e9
nh3_umol = ppb('NH3') * 1e-9 * rho / 28.97e-3 * 1e6; so4_umol = ug(SPEC['SO4']) / 96.06
FEAT = {**{f'f_{k}': fine[k] / mfine for k in SPEC}, 'log_mfine': np.log10(mfine + 0.05), 'log_mcoarse': np.log10(coarse + 0.5),
        'log_N': np.log10(ntot + 10), 'f_aitken': R[:, ix['num_a2_SRF']] * rho * 1e-6 / (ntot + 1e-3),
        'RH': R[:, ix['RHREFHT']], 'T': R[:, ix['TREFHT']], 'log_NH3': np.log10(ppb('NH3') + 0.01),
        'log_HNO3': np.log10(ppb('HNO3') + 0.001), 'log_SO2': np.log10(ppb('SO2') + 0.001),
        'log_NH3_2SO4': np.log10((nh3_umol + 1e-4) / (2 * so4_umol + 1e-4))}
X = np.stack([FEAT[n] for n in names], 1).astype(np.float32)
Zs = (X - mu) / sd * wfeat
day_of = np.repeat(np.arange(nd), len(cells)); cell_of = np.tile(cells, nd)
ok = np.isfinite(Zs).all(1)

# ---------------- selection -----------------------------------------------------------------------
ocean_like = fine['NCL'] / mfine > 0.5
sel, strata = [], []
tot_area = AREA.sum() * nd; tot_pop = pop.sum() * nd
cl = np.nan_to_num(R[:, ix['CLDLOW']]); tot_cld = (AREA[cell_of] * cl).sum() * ncell / len(cells)   # approx normaliser
for rname, (mask, metric, kind) in REGIONS.items():
    incell = mask[cell_of]
    for sname, mons in SEAS.items():
        m = ok & incell & np.isin(months[day_of], mons)
        if kind == 'land': m &= ~ocean_like
        if kind == 'ocean': m &= ocean_like | (fine['NCL'] / mfine > 0.2)
        idx = np.flatnonzero(m)
        if len(idx) < 20:
            print('skip', rname, sname, len(idx)); continue
        w_area = (AREA[cell_of[idx]]).sum() / tot_area; w_pop = pop[cell_of[idx]].sum() / tot_pop
        w_cld = (AREA[cell_of[idx]] * cl[idx]).sum()
        strata.append({'region': rname, 'season': sname, 'n_samples': int(len(idx)), 'area_time_frac': float(w_area),
                       'pop_time_frac': float(w_pop), 'lowcloud_area_time': float(w_cld)})
        val = (dust_tot if metric == 'dust' else mfine)[idx]
        for kind_s, sub in (('typical', idx), ('polluted', idx[val >= np.percentile(val, 90)])):
            med = np.median(Zs[sub], 0); k = sub[np.argmin(((Zs[sub] - med) ** 2).sum(1))]
            c = int(cell_of[k]); dd = int(day_of[k])
            sel.append({'region': rname, 'season': sname, 'kind': kind_s, 'date': dates[dd], 'cell': c,
                        'lat': float(LATf[c]), 'lon': float(LONf[c]), 'metric': metric,
                        'metric_value': float((dust_tot if metric == 'dust' else mfine)[k]),
                        'metric_p50': float(np.percentile(val, 50)), 'metric_p90': float(np.percentile(val, 90)),
                        'pm_fine': float(mfine[k]), 'dust': float(dust_tot[k]), 'raw': R[k].tolist()})
tot_lc = sum(x['lowcloud_area_time'] for x in strata)
for s in strata:
    s['lowcloud_frac'] = s['lowcloud_area_time'] / tot_lc   # share among strata
# population weight of the 30 area regimes
cnt = RG['counts'].astype(float); popw = (cnt * pop[None, :]).sum(1); popw /= popw.sum()
json.dump({'raw_names': RAW, 'samples': sel, 'strata': strata, 'regime_pop_weight': popw.tolist(),
           'regime_order': RG['order'].tolist(), 'pop_total': float(pop.sum())}, open(f'{OUT}/strata.json', 'w'))
for s in sel:
    print(s['region'], s['season'], s['kind'], s['date'], round(s['lat'], 1), round(s['lon'], 1), s['metric'], round(s['metric_value'], 2), round(s['metric_p90'], 2))
print('done', len(sel), 'samples', len(strata), 'strata')
