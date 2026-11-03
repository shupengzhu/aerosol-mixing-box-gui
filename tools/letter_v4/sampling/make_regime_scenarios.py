#!/usr/bin/env python3
"""Box scenarios from MAM4 regime medoids (regimes.npz from mam4_regimes.py).

Usage: python make_regime_scenarios.py regimes.npz <scenario_dir> [--init external|internal] [--ranks 0-29]

Design (forcing is sampled, not the internally mixed state):
  * met: daily-mean T, RH, PS of the medoid (grid cell, day); diurnal cycle T +/- 4 K over land
    (1 K where sea salt dominates), RH from constant specific humidity (cap 0.95, floor 0.15),
    J_NO2 = 9e-3 cos(SZA)^0.8 from latitude, day of year and local solar time; start 06:00, 24 h.
  * initial particles (init=external): each MAM4 mode is split into source populations with the
    mode's size distribution (number shared by volume): fresh carbonaceous (BC+POM), secondary
    (SO4 as NH4HSO4 + SOA), sea salt, dust; coarse-mode sulfate is put on the coarse dust and sea
    salt. init=internal keeps each MAM4 mode as one internally mixed population (sensitivity).
  * sources: dilution towards the initial state with tau = 3 d; primary particles (fresh BC/POM,
    dust, sea salt) emitted at M/tau_e (tau_e = 3 d fine, 1 d coarse); H2SO4 = k_OH,SO2 [OH] SO2,
    HNO3 = k_OH,NO2 [OH] NO2 with [OH] = 4e6 cm-3 x J/Jmax; SOA vapours (VBS bins as the MAM4
    particle SOA) at M_SOA/3 d, shaped by J_NO2.
  * MAM4 sigma_g (stock cam6_3 physprops): accumulation 1.6, Aitken 1.6, coarse 1.2, primary carbon 1.6.
"""
import json, math, sys, os
import numpy as np

src, outdir = sys.argv[1], sys.argv[2]
INIT = sys.argv[sys.argv.index('--init') + 1] if '--init' in sys.argv else 'external'
R0, R1 = (int(x) for x in (sys.argv[sys.argv.index('--ranks') + 1] if '--ranks' in sys.argv else '0-29').split('-'))
os.makedirs(outdir, exist_ok=True)
STRATA = src.endswith('.json')
if STRATA:                                   # stratified samples from mam4_strata.py
    SJ = json.load(open(src)); raw_names = SJ['raw_names']
    ITEMS = []
    for i, smp in enumerate(SJ['samples']):
        sid_ = f"strat_{smp['region']}_{smp['season']}_{smp['kind']}"
        ITEMS.append((i, sid_, np.array(smp['raw'], float), smp['lat'], smp['lon'] % 360, smp['date'],
                      f"{smp['region']} {smp['season']} {smp['kind']}: {smp['metric']} {smp['metric_value']:.2f} (P50 {smp['metric_p50']:.2f}, P90 {smp['metric_p90']:.2f})",
                      {'region': smp['region'], 'season': smp['season'], 'kind': smp['kind']}))
else:
    Z = np.load(src, allow_pickle=True)
    raw_names = list(Z['raw_names']); lat, lon = Z['lat'], Z['lon']; nlon = len(lon)
    ITEMS = []
    for rank, k in enumerate(Z['order']):
        if R0 <= rank <= R1:
            cell = int(Z['cand_cell'][k, 0])
            ITEMS.append((rank, f'regime_{rank:02d}', Z['cand_raw'][k, 0], float(lat[cell // nlon]), float(lon[cell % nlon]), str(Z['cand_day'][k, 0]),
                          f'cluster {int(k)}, area-time weight {float(Z["weight"][k]):.4f}, low-cloud weight {float(Z["weight_lowcloud"][k]):.4f}',
                          {'cluster': int(k), 'weight': float(Z['weight'][k]), 'weight_lowcloud': float(Z['weight_lowcloud'][k])}))
SIG = {'a1': 1.6, 'a2': 1.6, 'a3': 1.2, 'a4': 1.6}
RHO = {'so4': 1770., 'soa': 1000., 'pom': 1000., 'bc': 1700., 'dst': 2600., 'ncl': 1900.}
TAU_D, TAU_F, TAU_C = 3 * 86400., 3 * 86400., 86400.
OH_MAX, K_SO2, K_NO2 = 4e6, 1.0e-12, 1.1e-11          # cm-3, cm3 s-1
NA = 6.022e23


def val(raw, name):
    return float(raw[raw_names.index(name)])


def sza_cos(latd, doy, hour):
    dec = math.radians(23.44) * math.sin(2 * math.pi * (284 + doy) / 365)
    h = math.radians(15 * (hour - 12)); la = math.radians(latd)
    return math.sin(la) * math.sin(dec) + math.cos(la) * math.cos(dec) * math.cos(h)


def esat(T):
    return 611.2 * math.exp(17.67 * (T - 273.15) / (T - 29.65))


def dgn_from(N, V, sg):
    """number-median diameter (um) from number (m-3) and volume (m3 m-3) of a lognormal."""
    if N <= 0 or V <= 0:
        return None
    return (6 * V / (math.pi * N) * math.exp(-4.5 * math.log(sg) ** 2)) ** (1 / 3) * 1e6


def so4_split(m):
    return {'PSO4': m * 96.06 / 115.11, 'PNH4': m * 19.05 / 115.11}


def frac(d):
    tot = sum(d.values())
    return {k: round(v / tot, 6) for k, v in d.items() if v > 0} if tot > 0 else {}


rows = []
for rank, sid0, raw, la, lo, date, note, meta in ITEMS:
    T0, PS = val(raw, 'TREFHT'), val(raw, 'PS')
    RH0 = min(max(val(raw, 'RHREFHT') / 100.0, 0.15), 0.95)           # MAM4 RHREFHT is in percent
    rho = PS / (287.05 * T0)
    ug = lambda v: max(val(raw, v + '_SRF'), 0.0) * rho * 1e9          # kg/kg -> ug m-3
    ncm3 = lambda v: max(val(raw, v + '_SRF'), 0.0) * rho * 1e-6       # 1/kg -> cm-3
    ppb = lambda g: max(val(raw, g + '_SRF'), 0.0) * 1e9
    soa = {i: {m: ug(f'soa{i}_a{m}') for m in (1, 2)} for i in range(1, 6)}
    modes = {
        'a1': {'so4': ug('so4_a1'), 'soa': sum(soa[i][1] for i in soa), 'pom': ug('pom_a1'), 'bc': ug('bc_a1'), 'dst': ug('dst_a1'), 'ncl': ug('ncl_a1'), 'N': ncm3('num_a1')},
        'a2': {'so4': ug('so4_a2'), 'soa': sum(soa[i][2] for i in soa), 'dst': ug('dst_a2'), 'ncl': ug('ncl_a2'), 'N': ncm3('num_a2')},
        'a3': {'so4': ug('so4_a3'), 'dst': ug('dst_a3'), 'ncl': ug('ncl_a3'), 'N': ncm3('num_a3')},
        'a4': {'pom': ug('pom_a4'), 'bc': ug('bc_a4'), 'N': ncm3('num_a4')}}
    soa_bins = {f'PSOA{i}P': soa[i][1] + soa[i][2] for i in soa}
    marine = sum(modes[m].get('ncl', 0) for m in modes) > 0.5 * sum(v for m in modes for s, v in modes[m].items() if s != 'N')

    def species(comp, mode_soa):
        """MAM4 component masses (ug m-3) -> SCRAM species masses."""
        out = {}
        if comp.get('so4'): out.update({s: out.get(s, 0) + v for s, v in so4_split(comp['so4']).items()})
        if comp.get('soa'):
            tot = sum(mode_soa.values())
            for s, v in mode_soa.items():
                if v > 0: out[s] = out.get(s, 0) + comp['soa'] * v / tot
        if comp.get('pom'): out['PPOAlP'] = comp['pom']
        if comp.get('bc'): out['PBC'] = comp['bc']
        if comp.get('dst'): out['PMD'] = comp['dst']
        if comp.get('ncl'): out['PNA'] = 0.3866 * comp['ncl']; out['PHCL'] = 0.6134 * comp['ncl']
        return out

    pops = []
    for m, d in modes.items():
        sg = SIG[m]
        vol = {s: d[s] * 1e-9 / RHO[s] for s in d if s != 'N'}            # m3 m-3
        V = sum(vol.values())
        if d['N'] < 1e-3 or V <= 0:
            continue
        dg = dgn_from(d['N'] * 1e6, V, sg)
        msoa = {f'PSOA{i}P': soa[i][1 if m == 'a1' else 2] for i in soa} if m in ('a1', 'a2') else {}
        if INIT == 'internal':
            groups = {'internal': list(vol)}
        elif m == 'a3':
            cores = [s for s in ('dst', 'ncl') if vol.get(s, 0) > 0] or ['so4']
            groups = {s: [s] for s in cores}
        else:
            groups = {'carbonaceous': ['bc', 'pom'], 'secondary': ['so4', 'soa'], 'sea salt': ['ncl'], 'dust': ['dst']}
        for gname, members in groups.items():
            comp = {s: d[s] for s in members if d.get(s, 0) > 0}
            if INIT != 'internal' and m == 'a3' and d.get('so4', 0) > 0 and gname in ('dst', 'ncl'):
                core_tot = sum(d.get(s, 0) for s in groups)                 # coarse sulfate on the cores by mass
                comp['so4'] = d['so4'] * d[gname] / core_tot
            gv = sum(comp[s] * 1e-9 / RHO[s] for s in comp)
            if gv <= 0:
                continue
            n_g = d['N'] * gv / V
            if n_g < 1e-3:
                continue
            spec = species(comp, msoa)
            pops.append({'name': f'{m} {gname}', 'number_cm3': round(n_g, 4), 'dg_um': round(dg, 5), 'sigma_g': sg,
                         'mass_fraction': frac(spec)})
    # ---- emissions: primary particles at M / tau_e ---------------------------------------
    bcpom = sum(modes[m].get(s, 0) for m in ('a1', 'a4') for s in ('bc', 'pom'))
    emis = []
    def add_emis(name, mass_ugm3_s, comp, dg, sg, rho_p):
        if mass_ugm3_s <= 0: return
        v1 = math.pi / 6 * (dg * 1e-6) ** 3 * math.exp(4.5 * math.log(sg) ** 2)     # m3 per particle
        n = mass_ugm3_s * 1e-9 / rho_p / v1 / 1e6                                   # cm-3 s-1
        emis.append({'name': name, 'number_cm3_s': n, 'dg_um': dg, 'sigma_g': sg, 'mass_fraction': comp})
    if bcpom > 0:
        bc = sum(modes[m].get('bc', 0) for m in ('a1', 'a4')); pom = bcpom - bc
        add_emis('fresh carbonaceous', bcpom / TAU_F, frac({'PBC': bc, 'PPOAlP': pom}), 0.05, 1.7, (bc * 1700 + pom * 1000) / bcpom)
    for s, sp, rho_p in (('dst', {'PMD': 1.0}, 2600.), ('ncl', {'PNA': 0.3866, 'PHCL': 0.6134}, 1900.)):
        fine = modes['a1'].get(s, 0) + modes['a2'].get(s, 0); coarse = modes['a3'].get(s, 0)
        dg_c = dgn_from(modes['a3']['N'] * 1e6, sum(modes['a3'].get(x, 0) * 1e-9 / RHO[x] for x in ('so4', 'dst', 'ncl')), 1.2) or 2.0
        add_emis(f'{"dust" if s == "dst" else "sea salt"} fine', fine / TAU_F, sp, 0.3, 1.6, rho_p)
        add_emis(f'{"dust" if s == "dst" else "sea salt"} coarse', coarse / TAU_C, sp, round(min(dg_c, 4.0), 4), 1.2, rho_p)
    # ---- vapour sources (24 h mean scaled by J_NO2 shape) -------------------------------
    doy = (int(date[5:7]) - 1) * 30.4 + int(date[8:10])
    hours = np.arange(0, 25)                                               # t = 0 at 06:00
    mu = np.array([max(sza_cos(la, doy, (6 + h) % 24), 0.0) for h in hours])
    jno2 = 9e-3 * mu ** 0.8
    sunlit = jno2.max() > 0
    shape_mean = jno2[:24].mean() / jno2.max() if sunlit else 0.0         # 24-h mean of J/Jmax
    nair = PS / (1.380649e-23 * T0) * 1e-6                                 # molec cm-3
    # rates below are the PEAK (noon) values; the box scales gas sources by J_NO2/max(J_NO2)
    p_sulf = K_SO2 * OH_MAX * ppb('SO2') * 1e-9 * nair / NA * 98.08 * 1e6 * 1e6 if sunlit else 0.0   # ug m-3 s-1
    p_hno3 = K_NO2 * OH_MAX * ppb('NO2') * 1e-9 * nair / NA * 63.01 * 1e6 * 1e6 if sunlit else 0.0
    soa_tot = sum(soa_bins.values())
    gas_src = {'SULF': p_sulf, 'HNO3': p_hno3}
    if soa_tot > 0 and sunlit:                                            # 24-h mean production = M_SOA / tau
        gas_src.update({s.replace('PSOA', 'SOA'): v / TAU_F / shape_mean for s, v in soa_bins.items() if v > 0})
    gas_src = {s: float(f'{v:.4e}') for s, v in gas_src.items() if v > 0}
    gppb2ug = lambda g, mw: ppb(g) * 1e-9 * PS / (8.314 * T0) * mw * 1e6
    gas_init = {'NH3': gppb2ug('NH3', 17.03), 'HNO3': gppb2ug('HNO3', 63.01), 'SULF': gppb2ug('H2SO4', 98.08)}
    gas_init = {s: float(f'{v:.4e}') for s, v in gas_init.items() if v > 1e-6}
    # ---- meteorology ----------------------------------------------------------------------
    amp = 1.0 if marine else 4.0
    Tt = T0 - amp * np.cos(2 * np.pi * (hours - 9) / 24)                   # min at 06:00 (h=0), max at 15:00
    q = RH0 * esat(T0)
    RHt = np.clip(q / np.array([esat(t) for t in Tt]), 0.15, 0.95)
    sid = sid0 + ('' if INIT == 'external' else '_internal')
    case = {'case': {'name': sid, 'scenario_id': sid, 'category': 'regime',
                     'reference': f'MAM4 stock control sample {date} lat {la:.1f} lon {lo:.1f} ({note})'},
            'kernel': {'isorropia': 'v23'},
            'environment': {'temperature_K': float(Tt[0]), 'pressure_Pa': PS, 'relative_humidity': float(RHt[0]),
                            'time_series': {'time_s': (hours * 3600.0).tolist(), 'temperature_K': [round(float(x), 3) for x in Tt],
                                            'relative_humidity': [round(float(x), 4) for x in RHt], 'jno2_s1': [round(float(x), 7) for x in jno2]},
                            'dilution': {'rate_s1': 1 / TAU_D, 'background': 'initial'}},
            'time': {'duration_s': 86400.0, 'dt_s': 600.0, 'output_every_s': 3600.0},
            'mixing': {'preset': 'NC6'},
            'size_bins': {'mode': 'log', 'n': 16, 'dmin_um': 0.005, 'dmax_um': 10.0},
            'initial': {'particles': {'mode': 'lognormal', 'modes': pops}, 'gas_ugm3': gas_init},
            'emissions': {'particles': {'mode': 'lognormal', 'modes': emis}, 'gas_ugm3_s': gas_src, 'gas_time_profile': 'jno2' if sunlit else None}}
    scen = {'id': sid, 'category': 'regime', 'name': {'en': sid, 'zh': sid},
            'summary': {'en': f'Regime {rank} from MAM4 sampling', 'zh': f'MAM4 抽样区域状态 {rank}'}, 'tags': ['regime'],
            'notes': [note],
            'case': case}
    json.dump(scen, open(f'{outdir}/{sid}.json', 'w'), indent=1)
    rows.append({'id': sid, 'rank': rank, **meta, 'date': date, 'lat': la, 'lon': lo, 'marine': bool(marine),
                 'N_init_cm3': round(sum(p['number_cm3'] for p in pops), 1),
                 'PM_fine_ugm3': round(sum(modes[m].get(s, 0) for m in ('a1', 'a2', 'a4') for s in modes[m] if s != 'N'), 3),
                 'n_pops': len(pops), 'gas_src': gas_src, 'T': round(T0, 1), 'RH': round(RH0, 3)})
json.dump(rows, open(f'{outdir}/{"strata" if STRATA else "regime"}_index{"" if INIT == "external" else "_internal"}.json', 'w'), indent=1)
for r in rows:
    print(r['id'], r['date'], r['lat'], r['lon'], 'marine' if r['marine'] else 'land', r['N_init_cm3'], r['PM_fine_ugm3'], r['n_pops'], r['T'], r['RH'])
