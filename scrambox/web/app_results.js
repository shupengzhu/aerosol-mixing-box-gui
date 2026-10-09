// Runs, result views, experiment sets, start-up.
"use strict";
let exportPending = null;

// ================================================================== 5 run & results
function pageRun(page) {
  const nameInp = h("input", { type: "text", placeholder: t("run_name"), value: S.case ? S.case.case.name : "" });
  const pm = S.model === "partmc";
  const pmAvail = S.meta.partmc && Object.values(S.meta.partmc.available || {}).some(Boolean);
  S.pm = S.pm || Object.assign({}, (S.meta.partmc || {}).defaults || { n_part: 10000, del_t: 60, seed: 1, weighting_exponent: -1 });
  const pmField = (k, label, step) => field(label, h("input", { type: "number", step: step || "any", value: S.pm[k],
    onchange: ev => { S.pm[k] = parseFloat(ev.target.value); } }));
  const runBtn = h("button", { class: "btn primary", disabled: S.case && (!pm || pmAvail) ? null : true, onclick: async () => {
    try {
      postCfg();
      const body = { case: S.case, name: nameInp.value, model: S.model };
      if (pm) body.partmc = S.pm;
      const r = await api("runs", "POST", body);
      S.runId = r.id; S.results = null; toast(t("started") + ": " + r.id); await refreshRuns(); render();
    } catch (e) { toast(e.message, 6000); }
  } }, "▶ " + t(pm ? "run_partmc" : "run"));
  const lim = S.meta.limits;
  const pmCard = !pm ? null : h("div", { class: "card" }, h("h3", { style: { marginTop: 0 } }, t("partmc_opts")),
    pmAvail ? null : h("div", { class: "banner err" }, t("partmc_missing")),
    h("div", { class: "fields", style: { gridTemplateColumns: "1fr 1fr" } },
      pmField("n_part", t("n_part"), 1000), pmField("del_t", t("del_t"), 10), pmField("seed", t("seed"), 1), pmField("weighting_exponent", t("weighting"), 0.5)),
    h("p", { class: "small muted", style: { margin: "10px 0 4px" } }, t("partmc_hint")),
    h("p", { class: "small faint", style: { margin: 0 } }, t("partmc_not") +
      (lim ? ` ≤ ${lim.max_partmc_part} · ${S.lang === "zh" ? "粒子数×步数" : "particles × steps"} ≤ ${fmt(lim.max_partmc_cost)}` : "")));
  const left = h("div", {},
    h("div", { class: "card" }, h("div", { class: "field" }, nameInp), h("div", { class: "row", style: { marginTop: "8px" } }, runBtn)),
    pmCard,
    h("div", { class: "card" }, h("h3", { style: { marginTop: 0 } }, t("runs")), runList()));
  const right = h("div");
  page.append(h("div", { class: "side run-layout" }, left, right));
  if (S.runId) showRun(right);
  else right.append(h("div", { class: "card empty" }, t("no_runs"), h("div", { style: { marginTop: "12px" } }, exportMenu(null, false, t("export_select_run")))));
}

function runList() {
  if (!S.runs.length) return h("div", { class: "faint small" }, t("no_runs"));
  return h("div", { class: "runlist" }, S.runs.map(r => h("div", { class: "runitem" + (r.id === S.runId ? " on" : ""), onclick: () => { S.runId = r.id; S.results = null; S.tIndex = -1; render(); } },
    h("div", { class: "n" }, r.model === "partmc" ? h("span", { class: "chip accent", style: { marginRight: "6px", padding: "0 6px" } }, "PartMC") : null, r.name || r.id),
    h("div", { class: "m" }, h("span", { class: "dot " + r.state }), r.state, r.cpu_s ? `${fmt(r.cpu_s)} s` : "", h("span", { class: "spacer" }),
      h("button", { class: "btn ghost sm", title: t("delete"), onclick: async ev => { ev.stopPropagation(); await api("runs/" + r.id, "DELETE"); if (S.runId === r.id) S.runId = null; await refreshRuns(); render(); } }, "✕")),
    ["running", "queued"].includes(r.state) ? h("div", { class: "bar" }, h("i", { style: { width: ((r.progress || 0) * 100).toFixed(0) + "%" } })) : null)));
}

async function refreshRuns() { S.runs = await api("runs"); }

async function showRun(box) {
  const st = S.runs.find(r => r.id === S.runId) || {};
  if (st.state !== "done") {
    box.append(h("div", { class: "card" }, h("h2", {}, st.name || S.runId), h("div", { class: "row" }, h("span", { class: "dot " + st.state }), st.state, h("span", { class: "spacer" }), exportMenu(S.runId, false, t("export_wait"))),
      st.error ? h("pre", { class: "banner err mono", style: { whiteSpace: "pre-wrap" } }, st.error) : null,
      ["running", "queued"].includes(st.state) ? h("div", { class: "bar" }, h("i", { style: { width: ((st.progress || 0) * 100).toFixed(0) + "%" } })) : null));
    return;
  }
  if (!S.results || S.results.id !== S.runId) {
    box.append(h("div", { class: "card empty" }, "…"));
    S.results = await api(`runs/${S.runId}/results`);
    S.tIndex = S.results.time_s.length - 1;
    return render();
  }
  const R = S.results;
  const subs = [["overview", "overview"], ["size", "sizedist"], ["mix", "mixstate"], ["optics", "optics"], ["ccn", "ccn"], ["het", "hetchem"]];
  if (R.model === "partmc") subs.unshift(["particles", "particle_resolved"]);
  if (R.model !== "partmc" && S.resSub === "particles") S.resSub = "overview";
  box.append(h("div", { class: "card" },
    h("div", { class: "row" }, h("h2", {}, R.name), R.model === "partmc" ? h("span", { class: "chip accent" }, "PartMC" +
      (R.partmc_options ? ` · ${R.partmc_options.n_part} ${S.lang === "zh" ? "粒子" : "particles"}` : "")) : null, h("span", { class: "chip" }, R.kernel), h("span", { class: "chip" }, `${R.class_labels.length} ${t("classes_n")} × ${R.bounds_um.length - 1} ${t("bins_n")}`),
      h("span", { class: "chip" }, `cpu ${fmt(R.summary.cpu_s)} s`), h("span", { class: "spacer" }),
      exportMenu(R.id, true, null, R),
      h("button", { class: "btn sm", onclick: async () => { const c = await api(`runs/${R.id}/case`); await loadCase(c); toast(t("loaded")); } }, t("load_editor"))),
    R.notes && R.notes.length ? h("div", { class: "banner", style: { marginTop: "10px" } }, R.notes.join(" · ")) : null,
    h("div", { class: "subtabs", style: { marginTop: "12px" } }, subs.map(([id, k]) => h("button", { class: S.resSub === id ? "on" : "", onclick: () => { S.resSub = id; render(); } }, t(k))))));
  const body = h("div");
  box.append(body);
  ({ overview: resOverview, size: resSize, mix: resMix, optics: resOptics, ccn: resCcn, het: resHet, particles: resParticles })[S.resSub](body, R);
}

function exportMenu(runId, ready, reason, result) {
  const unavailable = exportPending ? t(exportPending === "save" ? "export_saving" : "export_packing") : reason;
  if (!ready || exportPending) return h("div", { class: "export-disabled" },
    h("button", { type: "button", class: "btn sm", disabled: true, title: unavailable }, t("export") + " ▾"),
    h("span", { class: "small muted" }, unavailable));
  const menu = h("details", { class: "export-menu" });
  const choices = [h("button", { type: "button", class: "btn ghost sm", onclick: () => exportRun(runId, "download", menu) }, t("export_download"))];
  if (S.mode === "trusted") choices.push(h("button", { type: "button", class: "btn ghost sm", onclick: () => exportRun(runId, "save", menu) }, t("export_save")));
  if (result) choices.push(h("button", { type: "button", class: "btn ghost sm", onclick: () => { menu.open = false; download(result, `${result.id}.json`); } }, t("export_json")));
  menu.append(h("summary", { class: "btn sm" }, t("export") + " ▾"), h("div", { class: "export-options" }, choices));
  menu.addEventListener("keydown", ev => { if (ev.key === "Escape") { menu.open = false; menu.querySelector("summary").focus(); } });
  menu.addEventListener("focusout", ev => { if (!menu.contains(ev.relatedTarget)) menu.open = false; });
  return menu;
}

async function exportRun(runId, destination, menu) {
  if (exportPending) return;
  exportPending = destination;
  menu.open = false;
  menu.setAttribute("aria-busy", "true");
  menu.querySelector("summary").textContent = t(destination === "save" ? "export_saving" : "export_packing");
  menu.querySelectorAll("button").forEach(button => { button.disabled = true; });
  try {
    if (destination === "save") {
      const response = await fetch("save-export", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ run_id: runId }) });
      if (response.redirected) throw new Error(t("export_session_expired"));
      const saveType = (response.headers.get("Content-Type") || "").toLowerCase();
      const result = saveType.includes("application/json") ? await response.json().catch(() => ({})) : {};
      if (!response.ok) throw new Error(result.error || (response.status === 401 || response.status === 403 ? t("export_session_expired") : t("export_failed")));
      if (!saveType.includes("application/json") || !result.ok) throw new Error(result.error || t("export_unexpected_response"));
      toast(t("export_saved") + ": " + (result.relative_path || result.filename), 6000);
      return;
    }
    const response = await fetch("api/runs/" + encodeURIComponent(runId) + "/export");
    if (response.redirected) throw new Error(t("export_session_expired"));
    const archiveType = (response.headers.get("Content-Type") || "").split(";")[0].trim().toLowerCase();
    if (!response.ok) {
      const error = archiveType.includes("json") ? await response.json().catch(() => ({})) : {};
      throw new Error(error.error || (response.status === 401 || response.status === 403 ? t("export_session_expired") : t("export_failed")));
    }
    if (archiveType !== "application/zip") throw new Error(t("export_unexpected_response"));
    const archive = await response.blob();
    if (!archive.size) throw new Error("Empty export");
    const disposition = response.headers.get("Content-Disposition") || "";
    const named = disposition.split("filename=")[1];
    const filename = (named ? named.split(";")[0].trim().replace(/^"|"$/g, "") : "scrambox-" + runId + ".zip").replace(/[\\/]/g, "_");
    const url = URL.createObjectURL(archive);
    const link = h("a", { href: url, download: filename });
    document.body.append(link); link.click(); link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    toast(t("export_downloaded") + ": " + filename, 5000);
  } catch (error) { toast(error.message, 6000); }
  finally {
    exportPending = null;
    if (menu.isConnected) {
      menu.removeAttribute("aria-busy");
      menu.querySelector("summary").textContent = t("export") + " ▾";
      menu.querySelectorAll("button").forEach(button => { button.disabled = false; });
    } else if (S.tab === "run") render();
  }
}

// ---------------------------------------------------------------- PartMC particle-resolved view
function resParticles(body, R) {
  const T = R.truth;
  const nOut = R.time_s.length;
  if (S.pIndex < 0 || S.pIndex >= nOut) S.pIndex = nOut - 1;
  const lab = h("span", { class: "mono" }, fmt(R.time_s[S.pIndex] / 3600) + " h");
  const sl = h("input", { type: "range", min: 0, max: nOut - 1, value: S.pIndex, style: { width: "260px" } });
  sl.oninput = () => { lab.textContent = fmt(R.time_s[+sl.value] / 3600) + " h"; };
  sl.onchange = () => { S.pIndex = +sl.value; S.particles = null; render(); };
  const info = h("span", { class: "muted small" });
  body.append(h("div", { class: "card" }, h("div", { class: "row" }, h("span", { class: "muted small" }, t("select_time")), sl, lab, h("span", { class: "spacer" }), info)));
  const charts = h("div", { class: "charts" });
  body.append(charts);
  if (T) {
    const th = T.time_s.map(x => x / 3600);
    // chi on the case grid, from the projected sections (same formula as tools/partmc_experiment.py)
    const iw = R.species.indexOf("PH2O");
    const chiSec = R.mass_ugm3.map(tb => {
      const secs = []; tb.forEach(row => row.forEach(cc => { const m = cc.filter((_, i) => i !== iw); const s = m.reduce((a, b) => a + b, 0); if (s > 0) secs.push(m); }));
      const tot = secs.reduce((a, m) => a + m.reduce((x, y) => x + y, 0), 0);
      let Da = 0; secs.forEach(m => { const s = m.reduce((x, y) => x + y, 0); const H = -m.reduce((a, v) => a + (v > 0 ? v / s * Math.log(v / s) : 0), 0); Da += s / tot * H; });
      const pa = secs[0] ? secs[0].map((_, j) => secs.reduce((a, m) => a + m[j], 0) / tot) : [];
      const Dg = Math.exp(-pa.reduce((a, v) => a + (v > 0 ? v * Math.log(v) : 0), 0));
      return Dg > 1 ? (Math.exp(Da) - 1) / (Dg - 1) : 1;
    });
    charts.append(chartCard(t("chi_title"), b => Charts.line(b, { xlabel: t("time_hours"), ymin: 0, ymax: 1, series: [
      { name: t("chi_particle"), x: th, y: T.chi, width: 2.6, color: "#1c2433" },
      { name: `${t("chi_sections")} (${R.class_labels.length} ${t("classes_n")})`, x: R.time_s.map(x => x / 3600), y: chiSec, dash: "6 4" }] })));
    const Ss = Object.keys(T.ccn_cm3);
    const post = R.post && R.post.ccn;
    charts.append(chartCard(t("ccn_truth") + ` (${fmt(R.time_s[nOut - 1] / 3600)} h)`, b => Charts.line(b, { logx: true, xlabel: "S (%)", ylabel: "cm⁻³", series: [
      { name: t("chi_particle"), x: Ss.map(Number), y: Ss.map(s => T.ccn_cm3[s][T.ccn_cm3[s].length - 1]), markers: true, width: 2.6, color: "#1c2433" },
      ...(post ? [{ name: t("chi_sections"), x: post.resolved[post.resolved.length - 1].spectrum.map(s => s.S_pct), y: post.resolved[post.resolved.length - 1].spectrum.map(s => s.ccn_cm3), dash: "6 4", markers: true }] : [])] })));
    if (T.bc_coating_ratio && T.bc_coating_ratio.some(v => v !== null)) {
      charts.append(chartCard(S.lang === "zh" ? "含 BC 颗粒：数浓度与包覆比（非 BC 干质量 / BC 质量）" : "BC-containing particles: number and coating ratio (non-BC dry mass / BC mass)",
        b => Charts.line(b, { xlabel: t("time_hours"), series: [{ name: "N_BC (cm⁻³)", x: th, y: T.n_bc_cm3 }, { name: S.lang === "zh" ? "包覆比" : "coating ratio", x: th, y: T.bc_coating_ratio }] })));
    }
  }
  const mapBC = h("div", { class: "card" }, h("h3", { style: { margin: 0 } }, t("bc_map")), h("div", { class: "empty" }, "…"));
  const mapK = h("div", { class: "card" }, h("h3", { style: { margin: 0 } }, t("kappa_map")), h("div", { class: "empty" }, "…"));
  charts.append(mapBC, mapK);
  const draw = P => {
    info.textContent = `${P.shown} / ${P.n_particles} ${t("particles_shown")}`;
    const opts = { x: P.d_um, w: P.n_cm3, xlabel: "D dry (µm)", wlabel: "N (cm⁻³)", nx: 40, ny: 25 };
    mapBC.lastChild.remove(); const b1 = h("div"); mapBC.append(b1);
    Charts.heat(b1, Object.assign({}, opts, { y: P.f_BC, ylabel: "BC mass fraction", ymin: 0, ymax: 1 }));
    mapK.lastChild.remove(); const b2 = h("div"); mapK.append(b2);
    Charts.heat(b2, Object.assign({}, opts, { y: P.kappa, ylabel: "κ", ymin: 0, ymax: Math.max(0.8, ...P.kappa.slice(0, 2000)) }));
  };
  if (S.particles && S.particles.id === R.id && S.particles.k === S.pIndex) draw(S.particles.data);
  else api(`runs/${R.id}/particles?k=${S.pIndex}`).then(P => { S.particles = { id: R.id, k: S.pIndex, data: P }; draw(P); })
    .catch(e => { mapBC.lastChild.textContent = e.message; });
}

function download(obj, name) {
  const a = h("a", { href: URL.createObjectURL(new Blob([JSON.stringify(obj)], { type: "application/json" })), download: name });
  document.body.append(a); a.click(); a.remove();
}

const th = R => R.time_s.map(x => x / 3600);
function chartCard(title, draw, extra) {
  const c = h("div", { class: "card" }, h("div", { class: "row" }, h("h3", { style: { margin: 0 } }, title), h("span", { class: "spacer" }), extra || null));
  const box = h("div"); c.append(box); setTimeout(() => draw(box), 0);
  return c;
}

function resOverview(body, R) {
  const iw = R.species.indexOf("PH2O"), nt = R.time_s.length, nc = R.class_labels.length;
  const N = R.number_m3.map(tb => tb.reduce((a, b) => a + b.reduce((x, y) => x + y, 0), 0) / 1e6);
  const Ncls = R.class_labels.map((_, c) => R.number_m3.map(tb => tb.reduce((a, b) => a + b[c], 0) / 1e6));
  const spMass = R.species.map((_, s) => R.mass_ugm3.map(tb => tb.reduce((a, b) => a + b.reduce((x, cc) => x + cc[s], 0), 0)));
  const act = R.species.map((s, i) => ({ s, i })).filter(o => o.i !== iw && Math.max(...spMass[o.i]) > 0);
  const gasIdx = R.gases.map((g, i) => ({ g, i })).filter(o => Math.max(...R.gas_ugm3.map(r => r[o.i])) > 1e-6);
  const gsel = h("select", {}, gasIdx.map(o => h("option", { value: o.i }, o.g)));
  const gasCard = chartCard(t("gases"), box => {
    const draw = () => {
      const chosen = [...gsel.selectedOptions].map(o => +o.value);
      Charts.line(box, { series: chosen.map(i => ({ name: R.gases[i], x: th(R), y: R.gas_ugm3.map(r => r[i]) })), xlabel: t("time_hours"), ylabel: "µg m⁻³" });
    };
    gsel.multiple = true; gsel.size = 1; gsel.style.height = "34px";
    [...gsel.options].slice(0, 4).forEach(o => o.selected = true);
    gsel.onchange = draw; draw();
  }, gsel);
  body.append(h("div", { class: "charts" },
    chartCard("N (cm⁻³)", box => Charts.line(box, { series: [{ name: "total", x: th(R), y: N, width: 2.5 }, ...Ncls.map((y, c) => ({ name: `C${c + 1}`, x: th(R), y, width: 1.4 }))], xlabel: t("time_hours"), ylabel: "cm⁻³" })),
    chartCard("PM dry (µg m⁻³)", box => Charts.line(box, { stacked: true, series: act.map(o => ({ name: o.s, x: th(R), y: spMass[o.i] })), xlabel: t("time_hours"), ylabel: "µg m⁻³" })),
    chartCard("PH2O (µg m⁻³)", box => Charts.line(box, { series: [{ name: "PH2O", x: th(R), y: spMass[iw], area: true }], xlabel: t("time_hours"), ylabel: "µg m⁻³" })),
    gasCard,
    chartCard("T / RH", box => Charts.line(box, { series: [{ name: "RH (%)", x: th(R), y: R.met.map(m => m[1] * 100) }, { name: "T − 273.15 (°C)", x: th(R), y: R.met.map(m => m[0] - 273.15) }], xlabel: t("time_hours") })),
    R.het ? chartCard(t("het_online") + " (µg m⁻³)", box => Charts.line(box, { xlabel: t("time_hours"), ylabel: "µg m⁻³", series: [
      { name: "N2O5 uptake (cum.)", x: th(R), y: R.het.map(r => r[3]) }, { name: "HNO3 produced (cum.)", x: th(R), y: R.het.map(r => r[4]) },
      { name: "ClNO2 produced (cum.)", x: th(R), y: R.het.map(r => r[5]) }] })) : null,
    R.het ? chartCard("k(N2O5) (s⁻¹) / φ(ClNO2)", box => Charts.line(box, { xlabel: t("time_hours"), series: [
      { name: "k (1/s)", x: th(R), y: R.het.map(r => r[0]) }, { name: "φ", x: th(R), y: R.het.map(r => r[2]) },
      { name: "[N2O5] (µg m⁻³)", x: th(R), y: R.het.map(r => r[6]) }] })) : null));
}

function timeSlider(R) {
  const lab = h("span", { class: "mono" }, fmt(R.time_s[S.tIndex] / 3600) + " h");
  const s = h("input", { type: "range", min: 0, max: R.time_s.length - 1, value: S.tIndex, style: { width: "260px" } });
  s.oninput = () => { lab.textContent = fmt(R.time_s[+s.value] / 3600) + " h"; };
  s.onchange = () => { S.tIndex = +s.value; render(); };
  return h("div", { class: "row" }, h("span", { class: "muted small" }, t("select_time")), s, lab);
}

function resSize(body, R) {
  const b = R.bounds_um, mid = b.slice(0, -1).map((v, i) => Math.sqrt(v * b[i + 1])), dl = b.slice(0, -1).map((v, i) => Math.log10(b[i + 1] / v));
  const ti = S.tIndex, iw = R.species.indexOf("PH2O");
  const dN = (k, c) => R.number_m3[k].map((row, bi) => (c === null ? row.reduce((a, x) => a + x, 0) : row[c]) / 1e6 / dl[bi]);
  const series = [{ name: "total", x: mid, y: dN(ti, null), width: 2.6, markers: true }, { name: "t = 0", x: mid, y: dN(0, null), dash: "5 4", color: "#8e99ad" },
    ...R.class_labels.map((l, c) => ({ name: `C${c + 1} ${l}`, x: mid, y: dN(ti, c), width: 1.3 }))];
  const act = R.species.map((s, i) => ({ s, i })).filter(o => o.i !== iw && R.mass_ugm3[ti].some(row => row.some(cc => cc[o.i] > 0)));
  body.append(h("div", { class: "card" }, timeSlider(R)), h("div", { class: "charts" },
    chartCard("dN/dlogD (cm⁻³)", box => Charts.line(box, { series, logx: true, logy: true, xlabel: "D (µm)", ylabel: "dN/dlogD" })),
    chartCard("dM/dlogD dry (µg m⁻³)", box => Charts.line(box, { stacked: true, logx: true, xlabel: "D (µm)", ylabel: "dM/dlogD",
      series: act.map(o => ({ name: o.s, x: mid, y: R.mass_ugm3[ti].map((row, bi) => row.reduce((a, cc) => a + cc[o.i], 0) / dl[bi]) })) }))));
}

function resMix(body, R) {
  const ti = S.tIndex, nc = R.class_labels.length, iw = R.species.indexOf("PH2O");
  const Ntot = R.number_m3.map(tb => tb.reduce((a, b) => a + b.reduce((x, y) => x + y, 0), 0));
  const frac = R.class_labels.map((_, c) => R.number_m3.map((tb, k) => tb.reduce((a, b) => a + b[c], 0) / (Ntot[k] || 1)));
  const act = R.species.map((s, i) => ({ s, i })).filter(o => o.i !== iw && R.mass_ugm3[ti].some(row => row.some(cc => cc[o.i] > 0)));
  const groups = R.class_labels.map((l, c) => `C${c + 1}`);
  body.append(h("div", { class: "card" }, timeSlider(R),
    h("div", { class: "chips", style: { marginTop: "8px" } }, R.class_labels.map((l, c) => h("span", { class: "chip" }, `C${c + 1} = ${l}`)))),
    h("div", { class: "charts" },
      chartCard(S.lang === "zh" ? "各类群数浓度占比" : "Number fraction per class", box => Charts.line(box, { stacked: true, ymax: 1, series: frac.map((y, c) => ({ name: `C${c + 1}`, x: th(R), y })), xlabel: t("time_hours") })),
      chartCard(S.lang === "zh" ? "各类群干质量组成 (µg m⁻³)" : "Dry mass by class and species (µg m⁻³)", box => Charts.bars(box, { stacked: true, groups,
        series: act.map(o => ({ name: o.s, values: R.class_labels.map((_, c) => R.mass_ugm3[ti].reduce((a, row) => a + row[c][o.i], 0)) })) }))));
}

function postMissing(body, R, key) {
  if (R.post && R.post[key]) return false;
  body.append(h("div", { class: "card" }, h("p", { class: "sub" }, S.lang === "zh" ? "该运行没有此后处理结果。可按当前“后处理”页设置计算：" : "This run has no such post-processing output. Compute it with the settings of the Post-processing tab:"),
    h("button", { class: "btn primary", onclick: async ev => {
      ev.target.disabled = true;
      try {
        const pp = S.case ? postCfg() : {};
        const cfg = JSON.parse(JSON.stringify(pp)); cfg[key] = Object.assign({}, cfg[key] || {}, { enabled: true });
        R.post = Object.assign(R.post || {}, await api(`runs/${R.id}/post`, "POST", cfg)); render();
      } catch (e) { toast(e.message, 5000); ev.target.disabled = false; }
    } }, t("compute_post"))));
  return true;
}

function resOptics(body, R) {
  if (postMissing(body, R, "optics")) return;
  const o = R.post.optics, tt = o.time_s.map(x => x / 3600);
  const last = o.resolved[o.resolved.length - 1], lastm = o.internal_mean[o.internal_mean.length - 1];
  body.append(h("div", { class: "charts" },
    chartCard("ext 550 nm (Mm⁻¹)", box => Charts.line(box, { series: [{ name: t("resolved"), x: tt, y: o.resolved.map(r => r.ext550_Mm) }, { name: t("internal_mean"), x: tt, y: o.internal_mean.map(r => r.ext550_Mm), dash: "6 4" },
      { name: "BC", x: tt, y: o.resolved.map(r => r.bc_ext550_Mm) }, { name: "dust", x: tt, y: o.resolved.map(r => r.dust_ext550_Mm) }], xlabel: t("time_hours"), ylabel: "Mm⁻¹" })),
    chartCard("SSA 550 nm", box => Charts.line(box, { series: [{ name: t("resolved"), x: tt, y: o.resolved.map(r => r.ssa550) }, { name: t("internal_mean"), x: tt, y: o.internal_mean.map(r => r.ssa550), dash: "6 4" }], xlabel: t("time_hours"), ymin: undefined })),
    chartCard(S.lang === "zh" ? "BC 吸收增强 E_abs" : "BC absorption enhancement E_abs", box => Charts.line(box, { series: [{ name: t("resolved"), x: tt, y: o.resolved.map(r => r.bc_abs_enhancement) }, { name: t("internal_mean"), x: tt, y: o.internal_mean.map(r => r.bc_abs_enhancement), dash: "6 4" }], xlabel: t("time_hours") })),
    chartCard(S.lang === "zh" ? "短波消光谱（终时刻）" : "SW extinction spectrum (last time)", box => Charts.line(box, { logx: true, xlabel: "λ (µm)", ylabel: "Mm⁻¹", series: [
      { name: t("resolved"), x: last.wavelength_sw_um, y: last.ext_sw_Mm, markers: true }, { name: t("internal_mean"), x: lastm.wavelength_sw_um, y: lastm.ext_sw_Mm, dash: "6 4" }] })),
    chartCard(S.lang === "zh" ? "各类群 550 nm 消光（终时刻）" : "550 nm extinction by class (last time)", box => Charts.bars(box, { groups: R.class_labels.map((_, c) => `C${c + 1}`), series: [{ name: "Mm⁻¹", values: last.ext550_by_class_Mm }] }))));
}

function resCcn(body, R) {
  if (postMissing(body, R, "ccn")) return;
  const c = R.post.ccn, tt = c.time_s.map(x => x / 3600);
  const last = c.resolved[c.resolved.length - 1], lastm = c.internal_mean[c.internal_mean.length - 1];
  body.append(h("div", { class: "charts" },
    chartCard(S.lang === "zh" ? "CCN 谱（终时刻）" : "CCN spectrum (last time)", box => Charts.line(box, { logx: true, xlabel: "S (%)", ylabel: "cm⁻³", series: [
      { name: t("resolved"), x: last.spectrum.map(s => s.S_pct), y: last.spectrum.map(s => s.ccn_cm3), markers: true },
      { name: t("internal_mean"), x: lastm.spectrum.map(s => s.S_pct), y: lastm.spectrum.map(s => s.ccn_cm3), dash: "6 4", markers: true }] })),
    chartCard("CCN(S) (cm⁻³)", box => Charts.line(box, { xlabel: t("time_hours"), ylabel: "cm⁻³", series: last.spectrum.map((s, k) => ({ name: `S=${s.S_pct}%`, x: tt, y: c.resolved.map(r => r.spectrum[k].ccn_cm3) })) })),
    chartCard(S.lang === "zh" ? "云滴数 Nd(w)（终时刻）" : "Droplet number Nd(w) (last time)", box => Charts.line(box, { xlabel: "w (m/s)", ylabel: "cm⁻³", series: [
      { name: t("resolved"), x: last.activation.map(a => a.w_m_s), y: last.activation.map(a => a.nd_cm3), markers: true },
      { name: t("internal_mean"), x: lastm.activation.map(a => a.w_m_s), y: lastm.activation.map(a => a.nd_cm3), dash: "6 4", markers: true }] })),
    chartCard(S.lang === "zh" ? "各类群活化比例（终时刻）" : "Activated fraction by class (last time)", box => Charts.bars(box, { groups: R.class_labels.map((_, k) => `C${k + 1}`),
      series: last.activation.map(a => ({ name: `w=${a.w_m_s}`, values: a.actfrac_by_class })) })),
    chartCard("kappa", box => Charts.bars(box, { groups: R.class_labels.map((_, k) => `C${k + 1}`), series: [{ name: "kappa", values: last.kappa_by_class.map(x => x || 0) }] }))));
}

function resHet(body, R) {
  if (postMissing(body, R, "hetchem")) return;
  const hc = R.post.hetchem, tt = hc.time_s.map(x => x / 3600), st = hc.states;
  const gases = Object.keys(st[0].rules);
  const charts = gases.map(g => {
    const rules = Object.keys(st[0].rules[g].resolved);
    return chartCard(`k(${g}) (s⁻¹)`, box => Charts.line(box, { logy: true, xlabel: t("time_hours"), ylabel: "s⁻¹", series: rules.flatMap((r, i) => [
      { name: `${r} · ${t("resolved")}`, x: tt, y: st.map(s => s.rules[g].resolved[r].k_s1), color: Charts.color(i) },
      { name: `${r} · ${t("internal_mean")}`, x: tt, y: st.map(s => s.rules[g].internal_mean[r].k_s1), color: Charts.color(i), dash: "6 4" }]) }));
  });
  const last = st[st.length - 1];
  const rows = [];
  gases.forEach(g => Object.keys(last.rules[g].resolved).forEach(r => {
    const a = last.rules[g].resolved[r], b = last.rules[g].internal_mean[r];
    rows.push(h("tr", {}, h("td", {}, `${g} · ${r}`), h("td", {}, fmt(a.k_s1)), h("td", {}, fmt(b.k_s1)), h("td", {}, fmt(a.gamma_eff)), h("td", {}, fmt(b.gamma_eff)),
      h("td", {}, b.k_s1 ? ((a.k_s1 / b.k_s1 - 1) * 100).toFixed(1) + "%" : "–"), h("td", {}, a.lifetime_h ? fmt(a.lifetime_h) : "–")));
  }));
  body.append(h("div", { class: "card" }, h("h3", { style: { marginTop: 0 } }, `SAD ${fmt(last.sad_um2_cm3)} µm² cm⁻³`),
    h("div", { class: "tablewrap" }, h("table", { class: "data" }, h("thead", {}, h("tr", {}, ["", "k " + t("resolved"), "k " + t("internal_mean"), "γeff " + t("resolved"), "γeff " + t("internal_mean"), "Δk", "τ (h)"].map(x => h("th", {}, x)))),
      h("tbody", {}, rows)))), h("div", { class: "charts" }, charts));
}

// ================================================================== 6 experiments
function pageExp(page) {
  const dims = S.meta.experiment_dimensions;
  const dimSel = h("select", {}, Object.keys(dims).map(k => h("option", { value: k }, k)), h("option", { value: "__custom" }, t("custom")));
  const custPath = h("input", { type: "text", placeholder: "path, e.g. environment.relative_humidity", style: { display: "none" } });
  const custVals = h("input", { type: "text", placeholder: "values, e.g. 0.5, 0.7, 0.9", style: { display: "none" } });
  dimSel.onchange = () => { const c = dimSel.value === "__custom"; custPath.style.display = custVals.style.display = c ? "" : "none"; };
  const nameInp = h("input", { type: "text", value: S.case ? S.case.case.name + " · " + (S.lang === "zh" ? "对比" : "comparison") : "" });
  const addBtn = h("button", { class: "btn", onclick: () => {
    let path = dimSel.value, vals;
    if (path === "__custom") { path = custPath.value.trim(); vals = custVals.value.split(",").map(s => s.trim()).filter(Boolean).map(s => { try { return JSON.parse(s); } catch { return s; } }); }
    else vals = dims[path];
    if (!path || !vals.length) return;
    vals.forEach(v => S.variants.push({ label: `${path.split(".").pop()}=${v}`, overrides: { [path]: v } }));
    render();
  } }, "+ " + t("add_variants"));
  const vrows = S.variants.map((v, i) => h("tr", {},
    h("td", {}, h("input", { type: "text", class: "txt", value: v.label, onchange: ev => v.label = ev.target.value })),
    h("td", {}, h("input", { type: "text", class: "txt", value: Object.entries(v.overrides).map(([k, x]) => `${k} = ${JSON.stringify(x)}`).join("; "),
      onchange: ev => { const o = {}; ev.target.value.split(";").forEach(p => { const [k, x] = p.split("="); if (k && x !== undefined) { try { o[k.trim()] = JSON.parse(x.trim()); } catch { o[k.trim()] = x.trim(); } } }); v.overrides = o; } })),
    h("td", { class: "rm" }, h("button", { class: "btn ghost sm", onclick: () => { S.variants.splice(i, 1); render(); } }, "✕"))));
  const runBtn = h("button", { class: "btn primary", disabled: S.case && S.variants.length ? null : true, onclick: async () => {
    try {
      const pp = postCfg();
      const r = await api("experiments", "POST", { name: nameInp.value, base: S.case, variants: S.variants, postprocess: pp });
      S.expId = r.id; S.expData = null; await refreshExps(); toast(t("started")); render();
    } catch (e) { toast(e.message, 5000); }
  } }, "▶ " + t("run_set"));
  const vsBtn = h("button", { class: "btn", onclick: () => {
    const pmo = S.pm || (S.meta.partmc || {}).defaults || {};
    S.variants = ["NC1", "NC6", "NC11"].map(nc => ({ label: `SCRAM ${nc}`, overrides: { "mixing.preset": nc } }))
      .concat([{ label: "PartMC", overrides: { model: "partmc", "partmc.n_part": pmo.n_part || 10000, "partmc.seed": pmo.seed || 1 } }]);
    render();
  } }, "⚖ " + t("vs_partmc"));
  page.append(h("div", { class: "card" }, h("h2", {}, t("experiments")), h("p", { class: "sub" }, t("exp_hint")),
    h("div", { class: "fields" }, field(t("name"), nameInp), field(t("dimension"), dimSel), field("", custPath), field("", custVals)),
    h("div", { class: "row", style: { marginTop: "10px" } }, addBtn, vsBtn, h("span", { class: "spacer" }), runBtn),
    S.variants.length ? h("div", { class: "tablewrap", style: { marginTop: "12px" } }, h("table", { class: "data edit" },
      h("thead", {}, h("tr", {}, h("th", {}, t("label")), h("th", {}, t("overrides")), h("th"))), h("tbody", {}, vrows))) : null));
  const list = h("div", { class: "runlist" }, S.exps.length ? S.exps.map(e => h("div", { class: "runitem" + (e.id === S.expId ? " on" : ""), onclick: () => { S.expId = e.id; S.expData = null; render(); } },
    h("div", { class: "n" }, e.name || e.id), h("div", { class: "m" }, h("span", { class: "dot " + e.state }), e.state, (e.variants || []).length + " " + t("variants"), h("span", { class: "spacer" }),
      h("button", { class: "btn ghost sm", onclick: async ev => { ev.stopPropagation(); await api("experiments/" + e.id, "DELETE"); if (S.expId === e.id) S.expId = null; await refreshExps(); render(); } }, "✕")))) : h("div", { class: "faint small" }, t("no_runs")));
  const right = h("div");
  page.append(h("div", { class: "side", style: { gridTemplateColumns: "280px minmax(0,1fr)" } }, h("div", { class: "card" }, h("h3", { style: { marginTop: 0 } }, t("experiments")), list), right));
  if (S.expId) showExp(right);
}

async function refreshExps() { S.exps = await api("experiments"); }

async function showExp(box) {
  if (!S.expData || S.expData.id !== S.expId) {
    box.append(h("div", { class: "card empty" }, "…"));
    S.expData = await api("experiments/" + S.expId);
    return render();
  }
  const E = S.expData;
  if (!E.variants) {
    box.append(h("div", { class: "card" }, h("h2", {}, E.status), h("ul", { class: "notes" }, (E.progress || []).map(v => h("li", {}, h("span", { class: "dot " + v.state }), " ", v.label, " — ", v.state)))));
    return;
  }
  const done = E.variants.filter(v => v.metrics);
  const metricDefs = [["N_cm3", "N (cm⁻³)", "time_h"], ["PM_dry", "PM dry (µg m⁻³)", "time_h"], ["PM1_dry", "PM1 dry (µg m⁻³)", "time_h"], ["water", "PH2O (µg m⁻³)", "time_h"]];
  if (done.length) {
    const m0 = done[0].metrics;
    Object.keys(m0.species).filter(s => s !== "PH2O").forEach(s => metricDefs.push(["species." + s, s + " (µg m⁻³)", "time_h"]));
    metricDefs.push(["chi_sections", "χ (" + t("chi_sections") + ")", "time_h"]);
    const pmRun = done.find(v => v.metrics.chi_particle);
    if (pmRun) {
      metricDefs.push(["chi_particle", "χ PartMC " + t("chi_particle"), "truth_time_h"]);
      Object.keys(pmRun.metrics.ccn_particle).forEach(s => metricDefs.push(["ccn_particle." + s, `CCN S=${s}% PartMC ${t("chi_particle")}`, "truth_time_h"]));
    }
    if (m0.ext550) metricDefs.push(["ext550", "ext 550 (Mm⁻¹)", "post_time_h"], ["ext550_internal", "ext 550 internal (Mm⁻¹)", "post_time_h"], ["ssa550", "SSA 550", "post_time_h"], ["bc_abs_enh", "BC E_abs", "post_time_h"]);
    if (m0.ccn) Object.keys(m0.ccn).forEach(s => metricDefs.push(["ccn." + s, `CCN S=${s}% (cm⁻³)`, "post_time_h"]));
    if (m0.nd) Object.keys(m0.nd).forEach(w => metricDefs.push(["nd." + w, `Nd w=${w} m/s (cm⁻³)`, "post_time_h"]));
    if (m0.k_het) Object.keys(m0.k_het).forEach(k => metricDefs.push(["k_het." + k, `k ${k} (s⁻¹)`, "post_time_h"]));
  }
  S.expMetric = S.expMetric && metricDefs.some(m => m[0] === S.expMetric) ? S.expMetric : metricDefs[0][0];
  const msel = h("select", { onchange: ev => { S.expMetric = ev.target.value; render(); } }, metricDefs.map(([k, l]) => h("option", { value: k, selected: k === S.expMetric ? true : null }, l)));
  const def = metricDefs.find(m => m[0] === S.expMetric);
  const getm = (v, k) => { const parts = k.split("."); return parts.length === 1 ? v.metrics[k] : (v.metrics[parts[0]] || {})[parts.slice(1).join(".")]; };
  const series = done.filter(v => v.metrics[def[2]] && getm(v, def[0]))
    .map(v => ({ name: v.label, x: v.metrics[def[2]], y: getm(v, def[0]) }));
  // grid metrics with a particle-resolved counterpart: overlay the PartMC truth as a dashed line
  const truthRun = done.find(v => v.metrics.chi_particle);
  if (truthRun && truthRun.metrics.truth_time_h) {
    let ty = null;
    if (def[0] === "chi_sections") ty = truthRun.metrics.chi_particle;
    else if (def[0].startsWith("ccn.")) {
      const sv = parseFloat(def[0].slice(4)), cp = truthRun.metrics.ccn_particle || {};
      const k = Object.keys(cp).find(k => Math.abs(parseFloat(k) - sv) < 1e-9);
      if (k) ty = cp[k];
    }
    if (ty) series.push({ name: truthRun.label + " · " + t("chi_particle"), x: truthRun.metrics.truth_time_h, y: ty, dash: "6 4" });
  }
  const finals = series.map(s => s.y[s.y.length - 1]);
  box.append(h("div", { class: "card" }, h("div", { class: "row" }, h("h2", {}, E.name), h("span", { class: "chip" }, E.status), h("span", { class: "spacer" }), field(t("metric"), msel)),
    E.variants.filter(v => v.error).map(v => h("div", { class: "banner err" }, `${v.label}: ${v.error}`))),
    chartCard(def[1], b => Charts.line(b, { series, xlabel: t("time_hours"), logy: def[0].startsWith("k_het") })),
    h("div", { class: "card" }, h("h3", { style: { marginTop: 0 } }, t("compare")), h("div", { class: "tablewrap" }, h("table", { class: "data" },
      h("thead", {}, h("tr", {}, h("th", {}, t("variants")), h("th", {}, t("final_value")), h("th", {}, t("vs_first")))),
      h("tbody", {}, series.map((s, i) => h("tr", {}, h("td", {}, s.name), h("td", {}, fmt(finals[i])), h("td", {}, finals[0] ? ((finals[i] / finals[0] - 1) * 100).toFixed(1) + "%" : "–"))))))));
}

// ================================================================== start-up and polling
async function poll() {
  try {
    const active = S.runs.some(r => ["running", "queued"].includes(r.state));
    const expActive = S.exps.some(e => ["running", "queued"].includes(e.state));
    if (active) { await refreshRuns(); if (S.tab === "run") { const cur = S.runs.find(r => r.id === S.runId); if (!cur || cur.state !== "done" || !S.results) render(); else renderRunListOnly(); } }
    if (expActive) { await refreshExps(); if (S.tab === "exp") { S.expData = null; render(); } }
  } catch (e) { /* server restarting */ }
  setTimeout(poll, 2000);
}
function renderRunListOnly() {
  const host = document.querySelector("#page-run .runlist");
  if (host) host.replaceWith(runList());
}

async function init() {
  document.documentElement.dataset.theme = (WEBSITE_EMBEDDED ? null : localStorage.getItem(storageKey("theme"))) || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
  $("#themeBtn").onclick = () => { const v = document.documentElement.dataset.theme === "dark" ? "light" : "dark"; document.documentElement.dataset.theme = v; localStorage.setItem(prefStorageKey("theme"), v); if (S.meta) render(); };
  document.querySelectorAll("#langSeg button").forEach(b => b.onclick = () => { S.lang = b.dataset.lang; localStorage.setItem(prefStorageKey("lang"), S.lang); document.documentElement.lang = S.lang === "zh" ? "zh-CN" : "en";
    if (S.meta) render(); else { document.querySelectorAll("[data-t]").forEach(e => e.textContent = t(e.dataset.t)); document.querySelectorAll("#langSeg button").forEach(x => x.classList.toggle("on", x.dataset.lang === S.lang)); showLogin(); } });
  document.querySelectorAll("#modelSeg button").forEach(b => b.onclick = () => {
    S.model = b.dataset.model;
    try { localStorage.setItem(prefStorageKey("model"), S.model); } catch (e) { /* storage unavailable */ }
    if (S.meta) { updateKernelChip(); render(); } else renderTabs();
  });
  const sess = await api("session");
  S.user = sess.user; S.mode = sess.mode;
  if (S.mode === "trusted") {
    S.lang = localStorage.getItem(prefStorageKey("lang")) || "zh";
    if (!PARAM_MODEL) S.model = localStorage.getItem(prefStorageKey("model")) || "scram";
    document.documentElement.lang = S.lang === "zh" ? "zh-CN" : "en";
    document.documentElement.dataset.theme = localStorage.getItem(prefStorageKey("theme")) || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
  }
  if (sess.mode === "shared" && !sess.user) { renderTabs(); $("#tabs").innerHTML = ""; return showLogin(); }
  userChip();
  S.meta = await api("meta");
  await reloadScenarios();
  await refreshRuns().catch(() => {});
  await refreshExps().catch(() => {});
  const first = S.scenarios.find(s => s.id === "urban_plume_aging") || S.scenarios[0];
  let draft = null;
  try { draft = JSON.parse(localStorage.getItem(caseStorageKey()) || "null"); } catch (e) { draft = null; }
  if (first) {
    S.scenario = await api("scenarios/" + first.id);
    S.describe = await api("describe", "POST", S.scenario.case).catch(e => ({ error: e.message }));
  }
  try { if (draft) await loadCase(draft); } catch (e) { draft = null; }
  if (!draft && S.scenario) await loadCase(S.scenario.case);
  render();
  poll();
}
setInterval(() => { if (S.dirty && S.case) { try { localStorage.setItem(caseStorageKey(), JSON.stringify(S.case)); S.dirty = false; } catch (e) { /* storage unavailable */ } } }, 3000);
init().catch(e => { document.querySelector("main").innerHTML = `<div class="banner err">${e.message}</div>`; });
