#!/usr/bin/env python3
"""Objective sampling of surface aerosol regimes from a stock CESM2/CAM6-MAM4 run.

Usage: python mam4_regimes.py <hist_dir> <case> <out_dir> [K=30] [day_stride=2] [years=2-5]

Reads daily h0 surface fields (bottom model layer), builds a mixing-state-relevant feature vector
for every (grid cell, day), clusters with area-weighted k-means (numpy only), and writes
  regimes.npz   centroids (raw units), weights, medoid candidates with full raw state, per-cell counts
  regimes.json  human-readable summary per cluster
Nothing in the case directory is written.
"""
import sys, os, glob, json, time
import numpy as np
from netCDF4 import Dataset

HIST, CASE, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
K = int(sys.argv[4]) if len(sys.argv) > 4 else 30
STRIDE = int(sys.argv[5]) if len(sys.argv) > 5 else 2
Y0, Y1 = (int(x) for x in (sys.argv[6] if len(sys.argv) > 6 else '2-5').split('-'))
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261011)

SPEC = {  # component: MAM4 tracers (a1 accumulation, a2 Aitken, a3 coarse, a4 primary carbon)
    'BC': ['bc_a1', 'bc_a4'], 'POM': ['pom_a1', 'pom_a4'],
    'SOA': [f'soa{i}_a{m}' for i in range(1, 6) for m in (1, 2)],
    'SO4': ['so4_a1', 'so4_a2', 'so4_a3'], 'DST': ['dst_a1', 'dst_a2', 'dst_a3'], 'NCL': ['ncl_a1', 'ncl_a2', 'ncl_a3']}
NUM = ['num_a1', 'num_a2', 'num_a3', 'num_a4']
GAS = ['SO2', 'H2SO4', 'HNO3', 'NH3', 'NO2', 'O3']
MET = ['TREFHT', 'RHREFHT', 'PS', 'PBLH', 'CLDLOW', 'U10', 'PM25_SRF', 'CCN3', 'CDNUMC', 'PRECT']
MVARS = sorted({v for vs in SPEC.values() for v in vs})
RAW = [v + '_SRF' for v in MVARS] + [n + '_SRF' for n in NUM] + [g + '_SRF' for g in GAS] + MET

files = []
for y in range(Y0, Y1 + 1):
    files += sorted(glob.glob(f'{HIST}/{CASE}.cam.h0.{y:04d}-*.nc'))
files = files[::STRIDE]
with Dataset(files[0]) as nc:
    lat, lon, gw = nc['lat'][:], nc['lon'][:], nc['gw'][:]
    missing = [v for v in RAW if v not in nc.variables]
assert not missing, missing
nlat, nlon, nd = len(lat), len(lon), len(files)
ncell = nlat * nlon
print(f'{nd} days x {ncell} cells, {len(RAW)} raw fields', flush=True)
raw = np.empty((nd, ncell, len(RAW)), np.float32)
dates = []
t0 = time.time()
for i, f in enumerate(files):
    with Dataset(f) as nc:
        for j, v in enumerate(RAW):
            x = nc[v]
            raw[i, :, j] = np.asarray(x[0, -1] if x.ndim == 4 else x[0], dtype=np.float32).ravel()   # 3-D: bottom layer
        dates.append(os.path.basename(f).split('.h0.')[1][:10])
    if i % 100 == 0:
        print(i, f'{time.time() - t0:.0f}s', flush=True)
ix = {v: j for j, v in enumerate(RAW)}
R = raw.reshape(nd * ncell, len(RAW))          # sample = day * ncell + cell

# ---- features -------------------------------------------------------------------------------
rho = R[:, ix['PS']] / (287.05 * R[:, ix['TREFHT']])                                      # kg m-3
ug = lambda vs: sum(R[:, ix[v + '_SRF']] for v in vs) * rho * 1e9                          # ug m-3
fine = {k: ug([v for v in vs if not v.endswith('_a3')]) for k, vs in SPEC.items()}
coarse = sum(ug([v for v in vs if v.endswith('_a3')]) for vs in SPEC.values())
mfine = sum(fine.values()) + 1e-6
ntot = sum(R[:, ix[n + '_SRF']] for n in NUM) * rho * 1e-6                                 # cm-3
ppb = lambda g: R[:, ix[g + '_SRF']] * 1e9
nh3_umol = ppb('NH3') * 1e-9 * rho / 28.97e-3 * 1e6                                         # umol m-3
so4_umol = ug(SPEC['SO4']) / 96.06
FEAT = {**{f'f_{k}': fine[k] / mfine for k in SPEC},
        'log_mfine': np.log10(mfine + 0.05), 'log_mcoarse': np.log10(coarse + 0.5),
        'log_N': np.log10(ntot + 10), 'f_aitken': R[:, ix['num_a2_SRF']] * rho * 1e-6 / (ntot + 1e-3),
        'RH': R[:, ix['RHREFHT']], 'T': R[:, ix['TREFHT']],
        'log_NH3': np.log10(ppb('NH3') + 0.01), 'log_HNO3': np.log10(ppb('HNO3') + 0.001),
        'log_SO2': np.log10(ppb('SO2') + 0.001), 'log_NH3_2SO4': np.log10((nh3_umol + 1e-4) / (2 * so4_umol + 1e-4))}
names = list(FEAT)
X = np.stack([FEAT[n] for n in names], 1).astype(np.float32)
del FEAT, fine
# composition fractions carry the mixing-state signal: weight them up relative to meteorology
wfeat = np.array([2.0 if n.startswith('f_') and n != 'f_aitken' else 1.0 for n in names], np.float32)
ok = np.isfinite(X).all(1)
mu, sd = X[ok].mean(0), X[ok].std(0) + 1e-6
Z = ((X - mu) / sd * wfeat).astype(np.float32)
wa = np.tile(np.repeat(gw, nlon), nd).astype(np.float64)                                    # area weight
print('finite', ok.mean(), flush=True)


def assign(Z, C, chunk=2_000_000):
    lab = np.empty(len(Z), np.int32); dmin = np.empty(len(Z), np.float32)
    cc = (C ** 2).sum(1)
    for s in range(0, len(Z), chunk):
        z = Z[s:s + chunk]
        d = cc[None, :] - 2 * z @ C.T
        lab[s:s + chunk] = d.argmin(1); dmin[s:s + chunk] = np.maximum(d.min(1) + (z ** 2).sum(1), 0)
    return lab, dmin


idx_ok = np.flatnonzero(ok)
Zo, wo = Z[idx_ok], wa[idx_ok] / wa[idx_ok].sum()
sub = rng.choice(len(Zo), 300_000, replace=False, p=wo)
C = [Zo[sub[0]]]
for _ in range(1, K):
    _, dm = assign(Zo[sub], np.array(C)); C.append(Zo[sub[rng.choice(len(sub), p=dm / dm.sum())]])
C = np.array(C, np.float32)
for it in range(60):
    lab, dmin = assign(Zo, C)
    W = np.bincount(lab, wo, K)
    Cn = np.stack([np.bincount(lab, wo * Zo[:, j], K) for j in range(Zo.shape[1])], 1) / np.maximum(W, 1e-30)[:, None]
    Cn[W == 0] = C[W == 0]
    shift = float(np.abs(Cn - C).max()); C = Cn.astype(np.float32)
    print('iter', it, f'shift {shift:.4f}', f'inertia {np.sum(dmin * wo):.4f}', flush=True)
    if shift < 2e-3:
        break
lab, dmin = assign(Zo, C)
weights = np.bincount(lab, wo, K)
cl = np.nan_to_num(R[idx_ok, ix['CLDLOW']])
wcld = np.bincount(lab, wo * cl, K); wcld /= wcld.sum()
order = np.argsort(-weights)
cand = np.full((K, 5), -1, np.int64)
for k in range(K):
    m = np.flatnonzero(lab == k)
    if len(m) >= 5:
        cand[k] = idx_ok[m[np.argsort(dmin[m])[:5]]]
Xo = X[idx_ok]
cent_raw = np.stack([np.bincount(lab, wo * Xo[:, j], K) for j in range(X.shape[1])], 1) / weights[:, None]
cnt = np.zeros((K, ncell), np.int32)
np.add.at(cnt, (lab, idx_ok % ncell), 1)
cand_raw = R[cand.ravel()].reshape(K, 5, len(RAW))
dates = np.array(dates)
np.savez_compressed(f'{OUT}/regimes.npz', names=names, raw_names=RAW, centroid=cent_raw, weight=weights, weight_lowcloud=wcld,
                    order=order, cand=cand, cand_raw=cand_raw, cand_day=dates[cand // ncell], cand_cell=cand % ncell,
                    counts=cnt, lat=lat, lon=lon, mu=mu, sd=sd, wfeat=wfeat, dates=dates)
summ = []
for r_, k in enumerate(order):
    c = cand[k, 0]; cell = c % ncell
    summ.append({'rank': r_, 'cluster': int(k), 'area_time_frac': round(float(weights[k]), 4), 'lowcloud_frac': round(float(wcld[k]), 4),
                 'medoid': {'date': str(dates[c // ncell]), 'lat': float(lat[cell // nlon]), 'lon': float(lon[cell % nlon])},
                 'centroid': {n: round(float(v), 3) for n, v in zip(names, cent_raw[k])}})
json.dump(summ, open(f'{OUT}/regimes.json', 'w'), indent=1)
print('done', flush=True)
