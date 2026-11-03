// Initial state, emissions and post-processing settings.
"use strict";

function drySpecies() { return S.grid.aerosol.filter(s => s !== "PH2O"); }
function binLabels() { const b = S.grid.bounds_um; return b.slice(0, -1).map((v, i) => `${fmt(v)}–${fmt(b[i + 1])}`); }
function classLabels() { return S.grid.class_labels.map((l, i) => `C${i + 1} ${l}`); }

// ================================================================== 3 inputs
function pageInputs(page) {
  if (!S.case) return page.append(h("div", { class: "card empty" }, t("scenario_hint")));
  const tabs = h("div", { class: "subtabs" }, [["initial", "initial"], ["emissions", "emissions"]].map(([id, k]) =>
    h("button", { class: S.inSub === id ? "on" : "", onclick: () => { S.inSub = id; render(); } }, t(k))));
  page.append(tabs);
  const sec = S.inSub;
  const rate = sec === "emissions";
  page.append(particleEditor(sec, rate));
  page.append(gasEditor(sec, rate));
  if (rate) page.append(profileEditor());
}

function emptyParticles(mode, rate) {
  const nb = S.grid.bounds_um.length - 1, nc = S.grid.class_labels.length;
  const kn = rate ? "number_cm3_s" : "number_cm3", km = rate ? "mass_ugm3_s" : "mass_ugm3";
  if (mode === "lognormal") return { mode, modes: [] };
  if (mode === "table") return { mode, [kn]: Array.from({ length: nc }, () => Array(nb).fill(0)), mass_fraction: Array.from({ length: nc }, () => ({})) };
  return { mode: "bin_mass", mixing: "internal", [km]: {}, [kn]: null };
}

function switchParticleMode(sec, next, rate) {
  const current = S.case[sec].particles;
  const mode = current.mode || "lognormal";
  if (mode === next) return;
  const drafts = S.particleModeDrafts[sec];
  drafts[mode] = clone(current);
  S.case[sec].particles = drafts[next] ? clone(drafts[next]) : emptyParticles(next, rate);
  S.dirty = true;
  S.caseDirty = true;
  render();
}

function particleEditor(sec, rate) {
  const p = S.case[sec].particles;
  const mode = p.mode || "lognormal";
  const kn = rate ? "number_cm3_s" : "number_cm3", km = rate ? "mass_ugm3_s" : "mass_ugm3";
  const seg = h("div", { class: "seg" }, [["lognormal", t("lognormal")], ["table", t("table_mode")], ["bin_mass", t("bin_mass")]].map(([v, l]) =>
    h("button", { class: mode === v ? "on" : "", onclick: () => switchParticleMode(sec, v, rate) }, l)));
  const card = h("div", { class: "card" }, h("div", { class: "row" }, h("h2", {}, t("particles")), h("span", { class: "spacer" }), h("span", { class: "muted small" }, t("particle_input")), seg));
  card.append(h("p", { class: "sub unit-note" }, t(rate ? "number_units_rate_hint" : "number_units_hint")));
  if (S.model === "partmc" && mode !== "lognormal") card.append(h("div", { class: "banner err", role: "alert" }, t("partmc_input_unsupported")));
  const sp = drySpecies();
  if (mode === "lognormal") {
    const modes = p.modes;
    const cols = [t("name"), t(rate ? "number_cm3_s" : "number_cm3"), t("dg"), t("sigma"), ...sp];
    const rows = modes.map((m, i) => h("tr", {},
      h("td", {}, h("input", { type: "text", class: "txt", value: m.name || "", onchange: ev => { m.name = ev.target.value; } })),
      ...[kn, "dg_um", "sigma_g"].map(k => h("td", {}, h("input", { type: "number", step: "any", value: m[k] ?? "", onchange: ev => { m[k] = parseFloat(ev.target.value); S.dirty = true; S.caseDirty = true; } }))),
      ...sp.map(s => {
        const v = (m.mass_fraction || {})[s] || 0;
        const inp = h("input", { type: "number", step: "any", min: 0, value: v || "" , placeholder: "0" });
        inp.addEventListener("change", () => { m.mass_fraction = m.mass_fraction || {}; const x = parseFloat(inp.value) || 0; if (x) m.mass_fraction[s] = x; else delete m.mass_fraction[s]; S.dirty = true; S.caseDirty = true; });
        return h("td", { class: v ? "" : "zero" }, inp);
      }),
      h("td", { class: "rm" }, h("button", { class: "btn ghost sm", onclick: () => { modes.splice(i, 1); render(); } }, "✕"))));
    card.append(h("p", { class: "sub" }, (S.lang === "zh" ? "每行一个对数正态模态；右侧各列为干质量分数（自动归一化），颗粒按组成自动归入类群。" :
      "One log-normal mode per row; the species columns are dry mass fractions (normalised), each mode is assigned to the class matching its composition.")),
      h("div", { class: "tablewrap" }, h("table", { class: "data edit" }, h("thead", {}, h("tr", {}, cols.map(c => h("th", {}, c)), h("th"))), h("tbody", {}, rows))),
      h("button", { class: "btn sm", style: { marginTop: "8px" }, onclick: () => { modes.push({ name: "mode " + (modes.length + 1), [kn]: rate ? 0.1 : 1000, dg_um: 0.1, sigma_g: 1.6, mass_fraction: { PSO4: 1 } }); render(); } }, "+ " + t("add_mode")));
  } else if (mode === "table") {
    const nb = S.grid.bounds_um.length - 1, nc = S.grid.class_labels.length;
    if (!p[kn] || p[kn].length !== nc || p[kn][0].length !== nb) { Object.assign(p, emptyParticles("table", rate)); }
    card.append(h("h3", {}, t("number_by_class") + ` (${t(rate ? "number_cm3_s" : "number_cm3")})`),
      editTable({ rows: classLabels(), cols: binLabels(), corner: "class \\ µm", get: (i, j) => p[kn][i][j], set: (i, j, v) => { p[kn][i][j] = v; } }),
      h("h3", {}, t("composition_by_class")),
      editTable({ rows: classLabels(), cols: sp, corner: "class", get: (i, j) => p.mass_fraction[i][sp[j]] || 0,
        set: (i, j, v) => { if (v) p.mass_fraction[i][sp[j]] = v; else delete p.mass_fraction[i][sp[j]]; } }));
  } else {
    const nb = S.grid.bounds_um.length - 1;
    const mass = p[km] = p[km] || {};
    Object.keys(mass).forEach(s => { if (!mass[s] || mass[s].length !== nb) mass[s] = Array(nb).fill(0); });
    card.append(h("div", { class: "fields", style: { marginBottom: "10px" } },
      field(t("mixing_of_input"), (() => {
        const s = h("select", { onchange: ev => { p.mixing = ev.target.value; S.dirty = true; S.caseDirty = true; render(); } },
          [["internal", t("internal")], ["by_group", t("by_group")], ["by_group_set", t("by_group_set")]].map(([v, l]) =>
            h("option", { value: v, selected: p.mixing === v ? true : null, disabled: v === "by_group_set" && !(p.sets && p.sets.length) ? true : null }, l)));
        return s;
      })())),
      p.mixing === "by_group_set" ? h("div", { class: "banner mixing-sets" }, t("by_group_set_hint"),
        h("div", { class: "mono" }, (p.sets || []).map((members, i) => `${i + 1}: ${members.join(" + ")}`).join(" | "))) : null,
      h("h3", {}, t("mass_by_bin") + (rate ? " s⁻¹" : "")),
      editTable({ rows: sp, cols: binLabels(), corner: "species \\ µm", get: (i, j) => (mass[sp[i]] || [])[j] || 0,
        set: (i, j, v) => { if (!mass[sp[i]]) mass[sp[i]] = Array(nb).fill(0); mass[sp[i]][j] = v; } }),
      h("h3", {}, (S.lang === "zh" ? "数浓度（可选，空=由质量与中值粒径推算）" : "Number (optional; empty = from mass at the bin mid diameter)")),
      editTable({ rows: [t(rate ? "number_cm3_s" : "number_cm3")], cols: binLabels(), corner: "",
        get: (_, j) => (p[kn] ? p[kn][j] : null), set: (_, j, v) => { if (!p[kn]) p[kn] = Array(nb).fill(0); p[kn][j] = v; } }));
  }
  return card;
}

function gasEditor(sec, rate) {
  const key = rate ? "gas_ugm3_s" : "gas_ugm3";
  const store = S.case[sec][key] = S.case[sec][key] || {};
  const card = h("div", { class: "card" });
  const filt = h("input", { type: "text", placeholder: t("filter"), style: { width: "200px" } });
  const nz = h("input", { type: "checkbox", checked: true });
  const holder = h("div");
  function draw() {
    const q = filt.value.trim().toLowerCase();
    const condensable = new Set(Object.values(S.meta.species[S.case.kernel.isorropia === "none" ? "noiso" : "iso"].gas_of));
    let gases = S.grid.gases.filter(g => !q || g.toLowerCase().includes(q));
    if (nz.checked && !q) gases = gases.filter(g => store[g] || condensable.has(g));
    gases.sort((a, b) => (condensable.has(b) - condensable.has(a)) || a.localeCompare(b));
    holder.innerHTML = "";
    holder.append(editTable({ rows: gases.map(g => condensable.has(g) ? g + " ●" : g), cols: [t(rate ? "gas_src" : "gas_conc")], corner: t("gases"),
      get: i => store[gases[i]] || 0, set: (i, _, v) => { if (v) store[gases[i]] = v; else delete store[gases[i]]; } }));
  }
  filt.addEventListener("input", draw); nz.addEventListener("change", draw);
  card.append(h("div", { class: "row" }, h("h2", {}, t("gases")), h("span", { class: "spacer" }), filt,
    h("label", { class: "check small" }, nz, t("nonzero_only") + " / ●"), ), holder);
  draw();
  return card;
}

function profileEditor() {
  const em = S.case.emissions;
  const card = h("div", { class: "card" }, h("h2", {}, t("emis_profile")));
  const prof = em.time_profile;
  const seg = h("div", { class: "seg" }, [["constant", t("constant")], ["custom", t("custom")]].map(([v, l]) => h("button", {
    class: (prof ? "custom" : "constant") === v ? "on" : "",
    onclick: () => { em.time_profile = v === "custom" ? { time_s: [0, S.case.time.duration_s], factor: [1, 1] } : null; render(); } }, l)));
  card.append(h("div", { class: "row" }, seg));
  const tab = (p, onRm) => h("div", { style: { marginTop: "10px", maxWidth: "380px" } }, editTable({ rows: p.time_s.map((_, i) => String(i + 1)), cols: [t("time_h"), t("factor")], corner: "#",
    get: (i, j) => j === 0 ? p.time_s[i] / 3600 : p.factor[i], set: (i, j, v) => { if (j === 0) p.time_s[i] = v * 3600; else p.factor[i] = v; },
    onRemoveRow: i => { if (p.time_s.length > 1) { p.time_s.splice(i, 1); p.factor.splice(i, 1); render(); } } }),
    h("button", { class: "btn sm", style: { marginTop: "6px" }, onclick: () => { p.time_s.push(p.time_s[p.time_s.length - 1] + 3600); p.factor.push(p.factor[p.factor.length - 1]); render(); } }, "+ " + t("add_row")));
  if (prof) card.append(tab(prof));
  const gp = em.gas_time_profile;
  const gmode = gp === "jno2" ? "jno2" : gp ? "custom" : "same";
  card.append(h("h3", {}, t("gas_profile")), h("div", { class: "seg" }, [["same", t("same_as_particles")], ["jno2", t("follow_jno2")], ["custom", t("custom")]].map(([v, l]) =>
    h("button", { class: gmode === v ? "on" : "", onclick: () => { em.gas_time_profile = v === "same" ? null : v === "jno2" ? "jno2" : { time_s: [0, S.case.time.duration_s], factor: [1, 1] }; render(); } }, l))));
  if (gp && typeof gp === "object") card.append(tab(gp));
  return card;
}

// ================================================================== 4 post-processing
function postCfg() {
  const pp = S.case.postprocess = S.case.postprocess || {};
  const pd = S.meta.post_defaults;
  pp.optics = Object.assign({ enabled: false, every: 1 }, pd.optics, pp.optics || {});
  pp.ccn = Object.assign({ enabled: false, every: 1 }, pd.ccn, pp.ccn || {});
  const hd = Object.assign({}, pd.hetchem); delete hd.species_gamma;
  pp.hetchem = Object.assign({ enabled: false, every: 1 }, hd, pp.hetchem || {});
  return pp;
}

function defaultGammaTable() {
  const tab = {}, chm = new Set(["PSO4", "PNH4", "PNO3", "PBC", "PSOA1P", "PSOA2P", "PSOA3P", "PSOA4P", "PSOA5P"]);
  for (const g in S.meta.post_defaults.gamma_cesm) {
    tab[g] = {};
    drySpecies().forEach(s => tab[g][s] = chm.has(s) ? S.meta.post_defaults.gamma_cesm[g] : 0);
  }
  return tab;
}

function pagePost(page) {
  if (!S.case) return page.append(h("div", { class: "card empty" }, t("scenario_hint")));
  const pp = postCfg();
  const o = pp.optics, c = pp.ccn, hc = pp.hetchem;
  page.append(h("div", { class: "card" }, h("div", { class: "row" }, h("h2", {}, t("optics")), h("span", { class: "spacer" }), toggle("postprocess.optics.enabled", t("enabled"))),
    h("p", { class: "sub" }, S.lang === "zh" ? "与 SCRAM CARMA_DiagnoseBulk 相同：RRTMG 14 个短波与 16 个长波波段，MAM4 kappa-Köhler 湿粒径，同一套 LUT。" :
      "Identical to SCRAM CARMA_DiagnoseBulk: RRTMG 14 SW and 16 LW bands, MAM4 kappa-Koehler wet radius, the same LUTs."),
    h("div", { class: "fields" },
      field(t("optics_mixing"), sel("postprocess.optics.mixing", [0, 1, 2, 3].map(v => ({ v, label: t("opt" + v) })))),
      field(t("dryshell"), sel("postprocess.optics.dryshell_selector", [{ v: 0, label: "0 sulfate-like (prod.)" }, { v: 1, label: "1 BrC-like" }, { v: 2, label: "2 org. fraction > 0.5" }])),
      field(t("every"), num("postprocess.optics.every", { step: 1 }))),
    h("div", { class: "row", style: { gap: "22px", marginTop: "12px" } }, toggle("postprocess.optics.perclass", t("perclass")), toggle("postprocess.optics.prognostic_moments", t("moments")))));
  page.append(h("div", { class: "card" }, h("div", { class: "row" }, h("h2", {}, t("ccn")), h("span", { class: "spacer" }), toggle("postprocess.ccn.enabled", t("enabled"))),
    h("p", { class: "sub" }, S.lang === "zh" ? "与 SCRAM scram_cloud_aerosol_export 相同：每个类群×粒径为一个模态，体积加权 kappa，线性质量倾斜的部分活化；给定过饱和度谱与 Abdul-Razzak & Ghan 上升气流活化。" :
      "Identical to SCRAM scram_cloud_aerosol_export: each class × bin is a mode with volume-weighted kappa and linear mass-tilt partial activation; fixed-supersaturation spectrum and Abdul-Razzak & Ghan activation for updrafts."),
    h("div", { class: "fields" }, field(t("supersats"), listInput("postprocess.ccn.supersaturations_pct")), field(t("updrafts"), listInput("postprocess.ccn.updrafts_m_s")),
      field(t("every"), num("postprocess.ccn.every", { step: 1 })))));
  const gt = hc.species_gamma = hc.species_gamma || defaultGammaTable();
  const gases = Object.keys(S.meta.post_defaults.gamma_cesm);
  const sp = drySpecies();
  sp.forEach(s => gases.forEach(g => { gt[g] = gt[g] || {}; if (gt[g][s] === undefined) gt[g][s] = 0; }));
  const gasChecks = h("div", { class: "row" }, gases.map(g => h("label", { class: "check" }, h("input", { type: "checkbox", checked: hc.gases.includes(g) ? true : null,
    onchange: ev => { hc.gases = ev.target.checked ? [...new Set([...hc.gases, g])] : hc.gases.filter(x => x !== g); } }), g)));
  const ruleChecks = h("div", { class: "row" }, [["cesm", "CESM (fixed γ)"], ["species", "species γ table"], ["bertram_thornton", "Bertram & Thornton (N2O5)"]].map(([r, l]) =>
    h("label", { class: "check" }, h("input", { type: "checkbox", checked: hc.rules.includes(r) ? true : null,
      onchange: ev => { hc.rules = ev.target.checked ? [...new Set([...hc.rules, r])] : hc.rules.filter(x => x !== r); } }), l)));
  page.append(h("div", { class: "card" }, h("div", { class: "row" }, h("h2", {}, t("hetchem")), h("span", { class: "spacer" }), toggle("postprocess.hetchem.enabled", t("enabled"))),
    h("p", { class: "sub" }, S.lang === "zh" ? "k = Σ sfc / (0.5 dm/Dg + 4/(c γ))，表面积与湿粒径同 SCRAM scram_hetchem_sad_export；同时给出类群分辨与内混参考。" :
      "k = Σ sfc / (0.5 dm/Dg + 4/(c γ)); surface and wet diameter as SCRAM scram_hetchem_sad_export; class-resolved and internal-mixture results side by side."),
    h("div", { class: "fields" }, field(t("het_gases"), gasChecks), field(t("het_rules"), ruleChecks)),
    h("div", { class: "row", style: { marginTop: "10px" } }, toggle("postprocess.hetchem.organic_coating", t("coating")), field(t("every"), num("postprocess.hetchem.every", { step: 1 }))),
    h("div", { class: "row" }, h("h3", {}, t("gamma_table")), h("span", { class: "spacer" }),
      h("button", { class: "btn sm", onclick: () => { hc.species_gamma = defaultGammaTable(); render(); } }, t("reset_cesm"))),
    editTable({ rows: sp, cols: gases, corner: "species \\ gas", get: (i, j) => gt[gases[j]][sp[i]], set: (i, j, v) => { gt[gases[j]][sp[i]] = v; } })));
}
