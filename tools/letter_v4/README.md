# Letter v4 tools (SCRAM vs PartMC error budget and configuration optimisation)

Sampling (run where the CESM2/MAM4 history and LandScan data live):
- sampling/mam4_regimes.py   area-weighted k-means background states (K = 30)
- sampling/mam4_strata.py    15 regions x 4 seasons x {typical, P90 polluted}, area / population / low-cloud weights
- sampling/make_regime_scenarios.py  box scenarios from the sampled forcing (--init external|internal)

Per scenario (tools/): partmc_experiment.py (PartMC + SCRAM NC1/6/11), stage2_regime.py (MAM4 projection,
parcel/ARG/FN chain, NC6 x8, NC9BC), wp2_optics.py, het_truth.py, aggregate_regimes.py.

Configuration optimisation (tools/): opt_configs.py (size grid x species groups x class bounds),
opt_screen.py (representation, PartMC projected; env OPT_SRC/OPT_OUT/OPT_SCRATCH),
opt_run.py (SCRAM simulation per configuration), opt_act_synth.py (732 synthetic activation cases).

Analysis and figures (letter_v4/): opt_merge.py + opt_analysis.py (tables), opt_stats.py (numbers in the text),
ensemble_figs.py, make_figs_v3.py (Fig. 2), make_figs_v4.py (Figs. 1 and 3), opt_figs.py (Fig. 4).
