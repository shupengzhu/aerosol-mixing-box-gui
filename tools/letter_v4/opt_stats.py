
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
def final_stats(an_dir, ens_csv, mam4_json, synth_json, e3c_json):
    L = pd.read_csv(f"{an_dir}/opt_long.csv"); C = pd.read_csv(f"{an_dir}/opt_configs.csv")
    E = pd.read_csv(ens_csv); M = json.load(open(mam4_json)); R = json.load(open(synth_json)); K = json.load(open(e3c_json))
    a = lambda x: (np.asarray(x, float) - 1) * 100
    def mp(x, sig=False):
        x = a(x); x = x[np.isfinite(x)]
        return (np.median(x), np.percentile(np.abs(x), 90)) if sig else (np.median(np.abs(x)), np.percentile(np.abs(x), 90))
    s = {}
    S = L[L.layer == "strata"]; G = L[L.layer == "regime"]
    b = S[S.cfg == "NC6@U1"]
    s["n_strat"] = b.id.nunique(); s["n_reg"] = G.id.nunique()
    s["R_med"], s["R_p90"] = mp(b["R0.1"])
    n1 = S[S.cfg == "NC1@U1"]; s["R1_med"], s["R1_p90"] = mp(n1["R0.1"])
    e = E[~E.internal.astype(bool) & (E.layer == "strata")]
    s["M_med"], s["M_p90"] = mp(e["M0.1"]); s["M_pos"] = int((e["M0.1"] > 1).sum()); s["M_n"] = int(e["M0.1"].notna().sum())
    pl = e[e.type == "polluted land"]; s["Mpl_med"], s["Mpl_p90"] = mp(pl["M0.1"])
    def cst(cfg, sub, col):
        r = C[(C.cfg == cfg) & (C.subset == sub)]
        return (float(r[col + "_med"].iloc[0]), float(r[col + "_p90"].iloc[0])) if len(r) and col + "_med" in r else (np.nan, np.nan)
    for cfg, k in (("NC6@U1", "b"), ("NC6@L2", "l2"), ("NC6@L4", "l4"), ("NC6@U8", "u8"), ("G2_3@U1", "g3"), ("G4_4@U1", "g4"), ("NC1@U1", "n1"), ("G1_6@U1", "g16")):
        for sub, sk in (("strata", "s"), ("strata_lowcloud", "c"), ("polluted_land", "p"), ("regimes_area", "r")):
            for col in ("S0.1", "Sall", "TSall", "Eall_mbn", "Sd0.5", "E_mbn0.5", "E_fn0.5", "E_arg0.5"):
                s[f"{k}_{sk}_{col}_med"], s[f"{k}_{sk}_{col}_p90"] = cst(cfg, sub, col)
        r = C[(C.cfg == cfg) & (C.subset == "strata")]
        s[f"{k}_cost"] = {"NC1@U1": 0.11, "NC6@U8": 87.9}.get(cfg, float(r.cost_rel.iloc[0]) if len(r) else np.nan)
        s[f"{k}_tr"] = int(r.tracers.iloc[0]) if len(r) else 0
    # time mean (NC6@U1, strata): representation from ensemble table
    s["TR_med"], s["TR_p90"] = np.median(e["TR0.1"]) * 100, np.percentile(e["TR0.1"], 90) * 100
    s["TS_med"], s["TS_p90"] = np.median(e["TS0.1"]) * 100, np.percentile(e["TS0.1"], 90) * 100
    # activation (NC6@U1 strata)
    for w in ("0.2", "0.5", "1.0"):
        D = b[f"Dpos{w}"].astype(str).str.lower() == "true"
        s[f"nDneg{w}"] = int((~D).sum())
        s[f"Aarg{w}"] = float(b[f"A_arg{w}"].median())
        s[f"Afn_pos{w}"], s[f"Afn_neg{w}"] = float(b[f"A_fn{w}"][D].median()), float(b[f"A_fn{w}"][~D].median())
        s[f"Ambn_pos{w}"], s[f"Ambn_neg{w}"] = float(b[f"A_mbn{w}"][D].median()), float(b[f"A_mbn{w}"][~D].median())
    # MAM4 + ARG end to end, w = 0.5
    ma = [((M.get(i) or {}).get("nd_arg") or {}).get("0.5", np.nan) / t for i, t in zip(b.id, b["nd_t0.5"])]
    s["MA_med_signed"] = float(np.nanmedian(a(ma))); s["MA_abs_med"] = float(np.nanmedian(np.abs(a(ma))))
    s["Mp_med_signed"] = float(np.nanmedian(a(e["A_mam0.5"])))
    # synthetic
    rr = [r for r in R if r["truth"] > 1]
    Dp = np.array([r["smax_FN"] ** 4 - 16 * K["A"] ** 2 * K["alpha"] * r["w"] / (9 * K["G"]) >= 0 for r in rr])
    f = lambda k, m: float(np.median([r[k] / r["truth"] for r, d in zip(rr, Dp) if d == m]))
    g = lambda k, m: float(np.median([r["smax_" + k] / r["smax_true"] for r, d in zip(rr, Dp) if d == m]))
    s.update(syn_n=len(rr), syn_neg=int((~Dp).sum()), syn_fn_neg=f("FN", False), syn_mbn_neg=f("MBN", False), syn_cnt_neg=f("count_at_true_smax", False),
             syn_fn_pos=f("FN", True), syn_mbn_pos=f("MBN", True), syn_sfn_neg=g("FN", False), syn_smbn_neg=g("MBN", False))
    # E_abs and k_N2O5 (strata, BC > 0.1 outside dust and ocean)
    for cfg, k in (("NC6@U1", "b"), ("G3_20@U1", "g320"), ("G4_4@U1", "g4")):
        q = S[(S.cfg == cfg) & (S.bc.astype(float) > 0.1) & (~S.type.isin(["dust", "ocean"]))]
        s[f"eabs_{k}_med"], s[f"eabs_{k}_p90"] = mp(q.eabs_sim / q.eabs_t)
        s[f"eabs_{k}_n"] = int((q.eabs_sim / q.eabs_t).notna().sum())
        q2 = S[S.cfg == cfg]; s[f"kbt_{k}_med"], s[f"kbt_{k}_p90"] = mp(q2.kbt_sim)
    # cross-validation of the ranking (season halves)
    S2 = S.copy(); S2["half"] = np.where(S2.season.isin(["DJF", "JJA"]), "A", "B"); S2["ea"] = (S2.Sall - 1).abs()
    pv = S2.groupby(["cfg", "half"]).ea.median().unstack().dropna()
    s["cv_rho"] = float(spearmanr(pv.A, pv.B).statistic); s["cv_n"] = len(pv)
    return s
