// SCRAM Box web interface (vanilla JS, no build step).
"use strict";

// ------------------------------------------------------------------ utilities
const S = {
  lang: localStorage.getItem("sb_lang") || "zh",
  meta: null, scenarios: [], scenario: null, describe: null,
  case: null, grid: null, tab: "scenario", setupSub: "kernel", inSub: "initial",
  runs: [], runId: null, results: null, resSub: "overview", tIndex: -1,
  exps: [], expId: null, expData: null, variants: [], dirty: false,
};
const $ = (sel, root = document) => root.querySelector(sel);
const t = k => (I18N[k] ? I18N[k][S.lang === "zh" ? 1 : 0] : k);
const L = obj => (obj && typeof obj === "object" ? (obj[S.lang] || obj.en || "") : (obj || ""));
const clone = o => JSON.parse(JSON.stringify(o));
const fmt = v => Charts.fmt(v);

function h(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  for (const k in attrs || {}) {
    const v = attrs[k];
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") e.className = v;
    else if (k === "style" && typeof v === "object") Object.assign(e.style, v);
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else if (k === "html") e.innerHTML = v;
    else e.setAttribute(k, v === true ? "" : v);
  }
  for (const c of kids.flat(Infinity)) if (c !== null && c !== undefined && c !== false) e.append(c.nodeType ? c : document.createTextNode(c));
  return e;
}

function getPath(o, p) { return p.split(".").reduce((a, k) => (a == null ? a : a[k]), o); }
function setPath(o, p, v) {
  const ks = p.split("."); let c = o;
  ks.slice(0, -1).forEach(k => { if (c[k] == null || typeof c[k] !== "object") c[k] = {}; c = c[k]; });
  c[ks[ks.length - 1]] = v;
}
function deepMerge(a, b) {
  const out = clone(a);
  for (const k in b || {}) {
    if (b[k] && typeof b[k] === "object" && !Array.isArray(b[k]) && out[k] && typeof out[k] === "object" && !Array.isArray(out[k])) out[k] = deepMerge(out[k], b[k]);
    else out[k] = clone(b[k]);
  }
  return out;
}
function toast(msg, ms = 2600) {
  const e = $("#toast"); e.textContent = msg; e.classList.add("show");
  clearTimeout(e._t); e._t = setTimeout(() => e.classList.remove("show"), ms);
}
async function api(path, method = "GET", body) {
  const r = await fetch("/api/" + path, { method, headers: { "Content-Type": "application/json" }, body: body ? JSON.stringify(body) : undefined });
  const j = await r.json().catch(() => ({ error: r.statusText }));
  if (r.status === 401 && j.login && path !== "login") { showLogin(); throw new Error(j.error); }
  if (!r.ok) throw new Error(j.error || r.statusText);
  return j;
}

// ------------------------------------------------------------------ sign-in (shared server)
function showLogin(msg) {
  const main = document.querySelector("main");
  $("#tabs").innerHTML = "";
  const code = h("input", { type: "password", placeholder: t("access_code"), autocomplete: "current-password" });
  const name = h("input", { type: "text", placeholder: t("your_name"), value: localStorage.getItem("sb_name") || "", autocomplete: "username" });
  const err = h("div", { class: "banner err", style: { display: msg ? "" : "none" } }, msg || "");
  const go = async () => {
    try {
      const r = await api("login", "POST", { code: code.value, name: name.value });
      localStorage.setItem("sb_name", name.value.trim());
      S.user = r.user; location.reload();
    } catch (e) { err.textContent = e.message; err.style.display = ""; }
  };
  [code, name].forEach(i => i.addEventListener("keydown", ev => { if (ev.key === "Enter") go(); }));
  main.innerHTML = "";
  main.append(h("div", { class: "card", style: { maxWidth: "440px", margin: "60px auto" } },
    h("h2", {}, t("sign_in")), h("p", { class: "sub" }, t("sign_in_hint")), err,
    h("div", { class: "fields", style: { gridTemplateColumns: "1fr" } }, field(t("access_code"), code), field(t("your_name"), name, t("name_hint"))),
    h("div", { class: "row", style: { marginTop: "14px", justifyContent: "flex-end" } }, h("button", { class: "btn primary", onclick: go }, t("sign_in")))));
}

function userChip() {
  const host = document.querySelector(".top-actions");
  const old = $("#userChip"); if (old) old.remove();
  if (!S.user) return;
  host.prepend(h("span", { id: "userChip", class: "chip", title: t("workspace_hint") }, "👤 " + S.user,
    h("button", { class: "btn ghost sm", style: { height: "20px", padding: "0 4px" }, onclick: async () => { await api("logout", "POST", {}); location.reload(); } }, t("sign_out"))));
}
const debounce = (fn, ms) => { let id; return (...a) => { clearTimeout(id); id = setTimeout(() => fn(...a), ms); }; };

// ------------------------------------------------------------------ form helpers
function field(label, input, hint) {
  return h("div", { class: "field" }, h("label", {}, label), input, hint ? h("div", { class: "hint" }, hint) : null);
}
function num(path, opts = {}) {
  const v = getPath(S.case, path);
  const scale = opts.scale || 1;
  const i = h("input", { type: "number", step: opts.step || "any", value: v === null || v === undefined ? "" : +(v / scale).toPrecision(10), placeholder: opts.placeholder || "" });
  i.addEventListener("change", () => {
    if (i.value === "" && opts.nullable) { setPath(S.case, path, null); changed(path); return; }
    const x = parseFloat(i.value);
    if (!isFinite(x)) { i.classList.add("bad"); return; }
    i.classList.remove("bad"); setPath(S.case, path, x * scale); changed(path);
  });
  return i;
}
function sel(path, options, onchange) {
  const v = getPath(S.case, path);
  const s = h("select", {}, options.map(o => h("option", { value: JSON.stringify(o.v), selected: JSON.stringify(o.v) === JSON.stringify(v) ? true : null }, o.label)));
  s.addEventListener("change", () => { setPath(S.case, path, JSON.parse(s.value)); changed(path); if (onchange) onchange(); });
  return s;
}
function toggle(path, label, onchange) {
  const i = h("input", { type: "checkbox", checked: getPath(S.case, path) ? true : null });
  i.addEventListener("change", () => { setPath(S.case, path, i.checked); changed(path); if (onchange) onchange(); });
  return h("label", { class: "check" }, h("span", { class: "switch" }, i, h("span")), label);
}
function listInput(path, opts = {}) {
  const v = getPath(S.case, path) || [];
  const i = h("input", { type: "text", value: v.join(", ") });
  i.addEventListener("change", () => {
    const xs = i.value.split(/[,\s]+/).filter(Boolean).map(Number);
    if (xs.some(x => !isFinite(x))) { i.classList.add("bad"); return; }
    i.classList.remove("bad"); setPath(S.case, path, xs); changed(path);
  });
  return i;
}

// Editable numeric table with pre-generated row and column names.
// get(r, c) -> number; set(r, c, value)
function editTable({ rows, cols, get, set, corner = "", rowExtra, onRemoveRow, colTitle }) {
  const thead = h("tr", {}, h("th", {}, corner), cols.map((c, j) => h("th", { title: colTitle ? colTitle(j) : null }, c)), onRemoveRow ? h("th") : null);
  const body = rows.map((r, i) => h("tr", {},
    h("td", {}, r),
    cols.map((_, j) => {
      const v = get(i, j);
      const inp = h("input", { type: "number", step: "any", value: v === null || v === undefined ? "" : +(+v).toPrecision(8) });
      const td = h("td", { class: +v === 0 ? "zero" : "" }, inp);
      inp.addEventListener("change", () => {
        const x = inp.value === "" ? 0 : parseFloat(inp.value);
        if (!isFinite(x)) { inp.classList.add("bad"); return; }
        inp.classList.remove("bad"); td.className = x === 0 ? "zero" : ""; set(i, j, x); S.dirty = true;
      });
      return td;
    }),
    onRemoveRow ? h("td", { class: "rm" }, h("button", { class: "btn ghost sm", onclick: () => onRemoveRow(i) }, "✕")) : null,
    rowExtra ? rowExtra(i) : null));
  return h("div", { class: "tablewrap" }, h("table", { class: "data edit" }, h("thead", {}, thead), h("tbody", {}, body)));
}

// ------------------------------------------------------------------ state changes
const GRID_KEYS = ["kernel", "species", "mixing", "size_bins", "processes.thermodynamics"];
const refreshGrid = debounce(async () => {
  try {
    S.grid = await api("grid", "POST", S.case);
    if (["setup", "inputs"].includes(S.tab)) render();
  } catch (e) { toast(e.message); }
}, 250);
function changed(path) {
  S.dirty = true;
  if (GRID_KEYS.some(k => path.startsWith(k))) refreshGrid();
  updateKernelChip();
}
function updateKernelChip() {
  if (!S.case) return;
  const k = S.case.kernel.isorropia;
  $("#kernelChip").textContent = `${k === "v23" ? "ISORROPIA II 2.3" : k === "v1" ? "ISORROPIA v1" : "no-ISO"} · ${S.case.mixing.preset} · ${S.case.case.name || ""}`;
}
async function loadCase(c) {
  S.case = deepMerge(S.meta.defaults, c);
  ["initial", "emissions"].forEach(sec => {
    if (!S.case[sec].particles) S.case[sec].particles = { mode: "lognormal", modes: [] };
  });
  S.grid = await api("grid", "POST", S.case);
  updateKernelChip();
}

// ------------------------------------------------------------------ tabs
const TABS = [["scenario", "tab_scenario"], ["setup", "tab_setup"], ["inputs", "tab_inputs"], ["post", "tab_post"], ["run", "tab_run"], ["exp", "tab_exp"]];
function renderTabs() {
  const nav = $("#tabs"); nav.innerHTML = "";
  TABS.forEach(([id, key], i) => nav.append(h("button", { class: "tab" + (S.tab === id ? " active" : ""), onclick: () => { S.tab = id; render(); } },
    h("span", { class: "num" }, String(i + 1)), t(key))));
  document.querySelectorAll("[data-t]").forEach(e => e.textContent = t(e.dataset.t));
  document.querySelectorAll("#langSeg button").forEach(b => b.classList.toggle("on", b.dataset.lang === S.lang));
}
function render() {
  renderTabs();
  document.querySelectorAll(".page").forEach(p => p.classList.remove("active"));
  const page = $("#page-" + S.tab); page.classList.add("active"); page.innerHTML = "";
  ({ scenario: pageScenario, setup: pageSetup, inputs: pageInputs, post: pagePost, run: pageRun, exp: pageExp })[S.tab](page);
}

// ================================================================== 1 scenario
function pageScenario(page) {
  const cats = S.meta.categories;
  const combo = h("div", { class: "combo" });
  const input = h("input", { type: "text", placeholder: t("search_scenario"), value: S.scenario ? L(S.scenario.name) : "" });
  const list = h("div", { class: "list" });
  combo.append(h("span", { class: "icon" }, "⌕"), input, list);
  let hl = 0, shown = [];
  const esc = s => s.replace(/[&<>]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" })[c]);
  const mark = (s, q) => !q ? esc(s) : esc(s).replace(new RegExp("(" + q.replace(/[.*+?^${}()|[\]\\]/g, "\\$&") + ")", "ig"), "<mark>$1</mark>");
  function fill() {
    const q = input.value.trim().toLowerCase();
    const terms = q.split(/\s+/).filter(Boolean);
    shown = S.scenarios.filter(s => {
      const hay = [s.id, s.name.en, s.name.zh, s.summary.en, s.summary.zh, (s.tags || []).join(" "), L(cats[s.category])].join(" ").toLowerCase();
      return terms.every(x => hay.includes(x));
    });
    list.innerHTML = ""; let cur = null;
    shown.forEach((s, i) => {
      if (s.category !== cur) { cur = s.category; list.append(h("div", { class: "grp" }, L(cats[cur]) || cur)); }
      const o = h("div", { class: "opt" + (i === hl ? " hl" : "") },
        h("div", { class: "t", html: mark(L(s.name), terms[0]) }), h("div", { class: "d", html: mark(L(s.summary), terms[0]) }));
      o.addEventListener("mousedown", ev => { ev.preventDefault(); choose(s.id); });
      list.append(o);
    });
    if (!shown.length) list.append(h("div", { class: "empty" }, "—"));
  }
  async function choose(id) {
    combo.classList.remove("open");
    S.scenario = await api("scenarios/" + id);
    S.describe = await api("describe", "POST", S.scenario.case).catch(e => ({ error: e.message }));
    render();
  }
  input.addEventListener("focus", () => { input.select(); hl = 0; fill(); combo.classList.add("open"); });
  input.addEventListener("input", () => { hl = 0; fill(); combo.classList.add("open"); });
  input.addEventListener("blur", () => setTimeout(() => combo.classList.remove("open"), 120));
  input.addEventListener("keydown", ev => {
    if (ev.key === "ArrowDown") { hl = Math.min(hl + 1, shown.length - 1); fill(); ev.preventDefault(); }
    else if (ev.key === "ArrowUp") { hl = Math.max(hl - 1, 0); fill(); ev.preventDefault(); }
    else if (ev.key === "Enter" && shown[hl]) { choose(shown[hl].id); input.blur(); }
    else if (ev.key === "Escape") input.blur();
  });

  page.append(h("div", { class: "card" }, h("h2", {}, t("scenario_library")), h("p", { class: "sub" }, t("scenario_hint")), combo));
  if (!S.scenario) return;
  const s = S.scenario, d = S.describe || {};
  const actions = h("div", { class: "row" },
    h("button", { class: "btn primary", onclick: async () => { await loadCase(s.case); toast(t("loaded")); S.tab = "setup"; render(); } }, t("load_editor")),
    h("button", { class: "btn", onclick: saveScenarioDialog }, t("save_my")),
    s.builtin ? null : h("button", { class: "btn danger", onclick: async () => { await api("scenarios/" + s.id, "DELETE"); S.scenario = null; await reloadScenarios(); render(); } }, t("delete")));
  page.append(h("div", { class: "card" },
    h("div", { class: "row" }, h("h2", {}, L(s.name)), h("span", { class: "chip accent" }, L(cats[s.category]) || s.category),
      (s.tags || []).map(x => h("span", { class: "chip" }, x)), h("span", { class: "spacer" }), actions),
    h("p", { style: { margin: "10px 0 6px" } }, L(s.summary)),
    s.case.case && s.case.case.reference ? h("div", { class: "small muted" }, t("reference") + ": " + s.case.case.reference) : null,
    s.notes && s.notes.length ? h("details", { style: { marginTop: "8px" } }, h("summary", {}, t("notes") + ` (${s.notes.length})`), h("ul", { class: "notes" }, s.notes.map(n => h("li", {}, n)))) : null));
  if (d.error) { page.append(h("div", { class: "banner err" }, d.error)); return; }
  page.append(describeCards(d));
}

function describeCards(d) {
  const e = d.environment, ip = d.initial_particles, pe = d.particle_emissions;
  const rng = (a, f = x => fmt(x)) => a[0] === a[1] ? f(a[0]) : `${f(a[0])} – ${f(a[1])}`;
  const stat = (k, v, unit) => h("div", { class: "stat" }, h("div", { class: "k" }, k), h("div", { class: "v" }, v, unit ? h("small", {}, " " + unit) : null));
  const kv = obj => {
    const ks = Object.keys(obj || {});
    if (!ks.length) return h("div", { class: "faint small" }, t("none"));
    return h("div", { class: "tablewrap", style: { maxHeight: "220px" } }, h("table", { class: "data" },
      h("tbody", {}, ks.sort((a, b) => obj[b] - obj[a]).map(k => h("tr", {}, h("td", {}, k), h("td", { class: "mono" }, fmt(obj[k])))))));
  };
  const proc = d.processes;
  return h("div", {},
    h("div", { class: "card" }, h("h3", { style: { marginTop: 0 } }, t("environment")),
      h("div", { class: "stats" },
        stat("T", rng(e.temperature_K), "K"), stat("RH", rng(e.relative_humidity, x => (100 * x).toFixed(0)), "%"),
        stat("P", rng(e.pressure_Pa, x => (x / 100).toFixed(0)), "hPa"), stat("J(NO2)", rng(e.jno2_s1), "s⁻¹"),
        stat(t("dilution"), rng(e.dilution_s1), ""), stat(t("duration"), fmt(d.time.duration_h), "h"),
        stat(t("substep"), fmt(d.time.substep_s), "s"), stat(t("kernel_calls"), d.time.n_kernel_calls, ""),
        stat("ISORROPIA", d.kernel === "v23" ? "II 2.3" : d.kernel === "v1" ? "v1" : "—", d.thermodynamics),
        stat(t("classes"), `${d.grid.n_class} × ${d.grid.n_size}`, d.grid.mixing),
        stat(t("coagulation") + " / " + t("condensation"), `${proc.coagulation ? "✓" : "✗"} / ${proc.condensation ? "✓" : "✗"}`, proc.nucleation ? "+ " + t("nucleation") : ""))),
    h("div", { class: "grid2" },
      h("div", { class: "card" }, h("h3", { style: { marginTop: 0 } }, t("background")),
        h("div", { class: "stats" }, stat("N", fmt(ip.number_cm3), "cm⁻³"), stat("PM1", fmt(ip.PM1), "µg m⁻³"),
          stat("PM2.5", fmt(ip.PM2_5), "µg m⁻³"), stat("PM10", fmt(ip.PM10), "µg m⁻³")),
        h("h3", {}, t("particles") + " (µg m⁻³)"), kv(ip.by_species), h("h3", {}, t("gases") + " (µg m⁻³)"), kv(d.initial_gas_ugm3)),
      h("div", { class: "card" }, h("h3", { style: { marginTop: 0 } }, t("emissions")),
        h("div", { class: "stats" }, stat("N", fmt(pe.number_cm3_per_s), "cm⁻³ s⁻¹"), stat("PM10", fmt(pe.PM10_per_s), "µg m⁻³ s⁻¹"),
          stat(t("emis_profile"), d.emission_profile, ""), stat(t("gas_profile"), d.gas_source_profile, "")),
        h("h3", {}, t("particles") + " (µg m⁻³ s⁻¹)"), kv(pe.by_species_per_s), h("h3", {}, t("gases") + " (µg m⁻³ s⁻¹)"), kv(d.gas_sources_ugm3_s))),
    d.notes && d.notes.length ? h("div", { class: "banner" }, d.notes.join(" · ")) : null);
}

async function reloadScenarios() { S.scenarios = await api("scenarios"); }

function saveScenarioDialog() {
  const dlg = $("#dlg"); dlg.innerHTML = "";
  const base = S.case ? S.case.case : {};
  const en = h("input", { type: "text", value: base.name || "" }), zh = h("input", { type: "text", value: "" });
  const sm = h("textarea", { rows: 3 }, base.description || "");
  dlg.append(h("h2", { style: { marginTop: 0 } }, t("save_my")),
    h("div", { class: "fields", style: { gridTemplateColumns: "1fr 1fr" } }, field(t("name_en"), en), field(t("name_zh"), zh)),
    h("div", { style: { marginTop: "10px" } }, field(t("summary"), sm)),
    h("div", { class: "row", style: { marginTop: "16px", justifyContent: "flex-end" } },
      h("button", { class: "btn", onclick: () => dlg.close() }, t("cancel")),
      h("button", { class: "btn primary", onclick: async () => {
        if (!S.case) await loadCase(S.scenario.case);
        const c = clone(S.case); c.case.name = en.value || zh.value || "my scenario";
        try {
          const s = await api("scenarios", "POST", { name: { en: en.value || zh.value, zh: zh.value || en.value }, summary: { en: sm.value, zh: sm.value },
            case: c, derived_from: S.scenario ? S.scenario.id : null });
          dlg.close(); await reloadScenarios(); S.scenario = await api("scenarios/" + s.id);
          S.describe = await api("describe", "POST", S.scenario.case); toast(t("saved")); render();
        } catch (e) { toast(e.message); }
      } }, t("save"))));
  dlg.showModal();
}

// ================================================================== 2 setup
function pageSetup(page) {
  if (!S.case) return page.append(h("div", { class: "card empty" }, t("scenario_hint")));
  const subs = [["kernel", "kernel_proc"], ["env", "env_time"], ["species", "species_groups"], ["classes", "classes"], ["bins", "size_bins"]];
  const nav = h("div", { class: "subnav" }, subs.map(([id, k]) => h("button", { class: S.setupSub === id ? "on" : "", onclick: () => { S.setupSub = id; render(); } }, t(k))));
  const body = h("div");
  ({ kernel: setupKernel, env: setupEnv, species: setupSpecies, classes: setupClasses, bins: setupBins })[S.setupSub](body);
  page.append(h("div", { class: "side" }, nav, body));
}

function setupKernel(body) {
  const kv = S.meta.kernels;
  const isoOpts = [{ v: "v23", label: t("iso_v23") + (kv.v23.available ? "" : " ✗") }, { v: "v1", label: t("iso_v1") + (kv.v1.available ? "" : " ✗") }, { v: "none", label: t("iso_none") }];
  const redis = S.meta.enums["processes.redistribution"].map(x => ({ v: x, label: x }));
  body.append(h("div", { class: "card" }, h("h2", {}, t("kernel_proc")),
    h("div", { class: "fields" },
      field(t("isorropia"), sel("kernel.isorropia", isoOpts, () => render())),
      field(t("thermo_scheme"), sel("processes.thermodynamics", [{ v: "dynamic", label: t("dynamic") }, { v: "equilibrium", label: t("equilibrium") }], () => render())),
      field(t("sulfate"), sel("processes.sulfate", [{ v: "dynamic", label: t("dynamic") }, { v: "equilibrium", label: t("equilibrium") }])),
      field(t("solver"), sel("processes.dynamic_solver", ["etr", "ros2", "euler"].map(x => ({ v: x, label: x.toUpperCase() })))),
      field(t("redistribution"), sel("processes.redistribution", redis)),
      field(t("wet_diam"), sel("processes.wet_diameter", [{ v: "isorropia", label: "ISORROPIA / closure" }, { v: "gerber", label: "Gerber" }])),
      field(t("density"), sel("processes.density.mode", [{ v: "variable", label: "variable" }, { v: "fixed", label: "fixed" }]))),
    h("hr", { class: "sep" }),
    h("div", { class: "row", style: { gap: "22px" } }, toggle("processes.coagulation", t("coagulation")), toggle("processes.condensation", t("condensation")),
      toggle("processes.nucleation.enabled", t("nucleation"), () => render()), toggle("processes.soa_photolysis", t("soa_photolysis"))),
    S.case.processes.nucleation.enabled ? h("div", { class: "fields", style: { marginTop: "12px" } },
      field(t("nucl_model"), sel("processes.nucleation.model", [{ v: "binary", label: "binary H2SO4-H2O" }, { v: "ternary", label: "ternary H2SO4-NH3-H2O" }]))) : null,
    h("details", { style: { marginTop: "14px" } }, h("summary", {}, t("advanced")),
      h("div", { class: "fields", style: { marginTop: "10px" } },
        ["bc_core_floor_um", "salt_core_floor_um", "dust_core_floor_um", "nucl_rate_cap_cm3_s", "nucl_h2so4_gas_fraction_max", "nucl_nh3_gas_fraction_max"]
          .map(k => field(k, num("advanced." + k))),
        field("fixed density (kg/m³)", num("processes.density.fixed_kgm3"))),
      h("div", { class: "row", style: { gap: "22px", marginTop: "10px" } },
        ["salt_rdb_bidirectional", "salt_fix_size_rdb", "dust_fix_size_rdb"].map(k => toggle("advanced." + k, k))))),
    hetCard(),
    h("div", { class: "card" }, h("h2", {}, "Kernel build"),
      ["v1", "v23"].map(v => h("details", {}, h("summary", {}, `${v}: ${kv[v].available ? "available" : "missing"}`), h("pre", { class: "mono small" }, kv[v].version)))));
}

function hetCard() {
  const on = getPath(S.case, "chemistry.n2o5_uptake.enabled");
  return h("div", { class: "card" },
    h("div", { class: "row" }, h("h2", {}, t("het_online")), h("span", { class: "spacer" }), toggle("chemistry.n2o5_uptake.enabled", t("enabled"), () => render())),
    h("p", { class: "sub" }, t("het_online_hint")),
    on ? h("div", { class: "fields" },
      field(t("n2o5_mode"), sel("chemistry.n2o5_uptake.mode", [{ v: "production", label: t("n2o5_prod") }, { v: "concentration", label: t("n2o5_conc") }], () => render())),
      getPath(S.case, "chemistry.n2o5_uptake.mode") === "concentration"
        ? field("N2O5 (µg m⁻³)", num("chemistry.n2o5_uptake.n2o5_ugm3"), t("n2o5_ts_hint"))
        : field("P(N2O5) (µg m⁻³ s⁻¹)", num("chemistry.n2o5_uptake.n2o5_production_ugm3_s"), "1 ppb/h ≈ 1.2e-3"),
      getPath(S.case, "chemistry.n2o5_uptake.mode") === "concentration" ? null : field("k_other (s⁻¹)", num("chemistry.n2o5_uptake.k_other_s1")),
      field(t("het_rules"), sel("chemistry.n2o5_uptake.rule", [{ v: "cesm", label: "CESM (γ on f_chm surface)" }, { v: "species", label: t("gamma_table") }])),
      field("γ (CESM)", num("chemistry.n2o5_uptake.gamma_cesm")),
      field("ClNO2", sel("chemistry.n2o5_uptake.clno2_yield", [{ v: "bertram_thornton", label: "φ = 1/(1+[H2O]/(483[Cl⁻]))" }, { v: "none", label: t("none") }]))) : null);
}

function setupEnv(body) {
  const c = S.case;
  const ts = c.environment.time_series;
  const defSub = c.kernel.isorropia !== "none" && c.processes.thermodynamics === "dynamic" ? 10 : 600;
  const card1 = h("div", { class: "card" }, h("h2", {}, t("environment")),
    h("div", { class: "fields" },
      field(t("temperature"), num("environment.temperature_K")), field(t("rh"), num("environment.relative_humidity")),
      field(t("pressure"), num("environment.pressure_Pa")), field(t("jno2"), num("environment.jno2_s1")),
      field(t("dilution"), num("environment.dilution.rate_s1"), "→ " + (S.lang === "zh" ? "背景 = 初始状态" : "background = initial state"))),
    h("hr", { class: "sep" }),
    h("label", { class: "check" }, h("span", { class: "switch" }, h("input", { type: "checkbox", checked: ts ? true : null, onchange: ev => {
      if (ev.target.checked) {
        const e = c.environment;
        c.environment.time_series = { time_s: [0, c.time.duration_s], temperature_K: [e.temperature_K, e.temperature_K], relative_humidity: [e.relative_humidity, e.relative_humidity],
          pressure_Pa: [e.pressure_Pa, e.pressure_Pa], jno2_s1: [e.jno2_s1 || 0, e.jno2_s1 || 0], dilution_s1: [e.dilution.rate_s1 || 0, e.dilution.rate_s1 || 0] };
      } else c.environment.time_series = null;
      changed("environment"); render();
    } }), h("span")), t("time_varying")));
  if (ts) {
    const cols = [["time_s", t("time_h"), 3600], ["temperature_K", "T (K)", 1], ["relative_humidity", "RH", 1], ["pressure_Pa", "P (Pa)", 1], ["jno2_s1", "J(NO2)", 1], ["dilution_s1", t("dilution"), 1], ["n2o5_production_ugm3_s", "P(N2O5)", 1], ["n2o5_ugm3", "N2O5 (µg m⁻³)", 1]];
    cols.forEach(([k]) => { if (!ts[k]) ts[k] = ts.time_s.map(() => (k === "pressure_Pa" ? c.environment.pressure_Pa : 0)); });
    card1.append(h("div", { style: { marginTop: "12px" } }, editTable({
      rows: ts.time_s.map((_, i) => String(i + 1)), cols: cols.map(x => x[1]), corner: "#",
      get: (i, j) => ts[cols[j][0]][i] / cols[j][2], set: (i, j, v) => { ts[cols[j][0]][i] = v * cols[j][2]; },
      onRemoveRow: i => { if (ts.time_s.length > 1) { cols.forEach(([k]) => ts[k].splice(i, 1)); render(); } } })),
      h("div", { class: "row", style: { marginTop: "8px" } }, h("button", { class: "btn sm", onclick: () => { cols.forEach(([k]) => ts[k].push(ts[k][ts[k].length - 1] + (k === "time_s" ? 3600 : 0))); render(); } }, "+ " + t("add_row"))));
  }
  const card2 = h("div", { class: "card" }, h("h2", {}, "Time"),
    h("div", { class: "fields" },
      field(t("duration"), num("time.duration_s", { scale: 3600 })), field(t("dt"), num("time.dt_s")),
      field(t("substep"), num("time.substep_s", { nullable: true, placeholder: `auto: ${defSub}` }), t("substep_hint")),
      field(t("output_every"), num("time.output_every_s"))));
  body.append(card1, card2);
}

function setupSpecies(body) {
  const lib = S.meta.species[S.case.kernel.isorropia === "none" ? "noiso" : "iso"];
  const g = S.grid;
  body.append(h("div", { class: "card" }, h("h2", {}, t("library_species")),
    h("div", { class: "chips", style: { marginBottom: "10px" } }, g.groups.map((x, i) => h("span", { class: "chip accent" }, `${i + 1}. ${x}`))),
    h("div", { class: "tablewrap" }, h("table", { class: "data" },
      h("thead", {}, h("tr", {}, [t("name"), t("group"), t("density_gcm3"), t("mw"), t("kappa"), t("cond_gas")].map(x => h("th", {}, x)))),
      h("tbody", {}, lib.aerosol.map(s => h("tr", {}, h("td", {}, s), h("td", {}, lib.group_of[s] || "—"), h("td", {}, fmt(lib.density_gcm3[s])),
        h("td", {}, fmt(lib.mw_gmol[s])), h("td", {}, fmt(lib.kappa[s])), h("td", {}, lib.gas_of[s] || "—"))))))));
  const custom = S.case.species.custom;
  const rows = custom.map((c, i) => h("tr", {},
    h("td", {}, h("input", { type: "text", class: "txt", value: c.name, onchange: ev => { c.name = ev.target.value.trim(); changed("species"); } })),
    h("td", {}, h("select", { onchange: ev => { c.group = ev.target.value; changed("species"); } }, g.groups.map(x => h("option", { value: x, selected: x === c.group ? true : null }, x)))),
    ...["density_gcm3", "mw_gmol", "kappa"].map(k => h("td", {}, h("input", { type: "number", step: "any", value: c[k] ?? "", onchange: ev => { c[k] = parseFloat(ev.target.value); changed("species"); } }))),
    h("td", { class: "rm" }, h("button", { class: "btn ghost sm", onclick: () => { custom.splice(i, 1); changed("species"); render(); } }, "✕"))));
  body.append(h("div", { class: "card" }, h("h2", {}, t("custom_species")), h("p", { class: "sub" }, t("custom_hint")),
    custom.length ? h("div", { class: "tablewrap" }, h("table", { class: "data edit" },
      h("thead", {}, h("tr", {}, [t("name"), t("group"), t("density_gcm3"), t("mw"), t("kappa"), ""].map(x => h("th", {}, x)))), h("tbody", {}, rows))) : null,
    h("div", { class: "row", style: { marginTop: "10px" } }, h("button", { class: "btn sm", onclick: () => {
      custom.push({ name: "X" + (custom.length + 1), group: g.groups[g.groups.length - 1], density_gcm3: 1.5, mw_gmol: 150, kappa: 0.1 });
      changed("species"); render();
    } }, "+ " + t("add_species")))));
}

function setupClasses(body) {
  const m = S.case.mixing, g = S.grid;
  const presets = ["NC1", "NC6", "NC11", "custom"].map(x => ({ v: x, label: x === "custom" ? t("custom_preset") : x }));
  const card = h("div", { class: "card" }, h("h2", {}, t("classes")),
    h("div", { class: "fields" }, field(t("mixing_preset"), sel("mixing.preset", presets, () => render()))));
  if (m.preset === "custom") card.append(h("div", { class: "fields", style: { marginTop: "12px" } },
    h("div", { class: "field" }, h("label", {}, t("external")), toggle("mixing.external", "")),
    field(t("primary_bounds"), listInput("mixing.primary_fraction_bounds")),
    field(t("secondary_bounds"), listInput("mixing.secondary_fraction_bounds"))));
  card.append(h("h3", {}, `${t("generated_classes")}: ${g.class_labels.length}`),
    h("div", { class: "classbox" }, g.classes.map((cls, i) => h("div", { class: "cls" }, h("div", { class: "h" }, `C${i + 1} · ${g.class_labels[i]}`),
      h("div", { class: "b" }, cls.map((b, j) => h("span", { class: b[0] > 0 || b[1] < 1 ? "hi" : "" }, `${g.groups[j]} ${b[0]}–${b[1]}`)))))));
  body.append(card);
}

function setupBins(body) {
  const sb = S.case.size_bins;
  const mode = sb.mode === "explicit" ? "explicit" : "log";
  const seg = h("div", { class: "seg" }, [["log", t("log_auto")], ["explicit", t("manual")]].map(([v, l]) => h("button", { class: mode === v ? "on" : "", onclick: () => {
    if (v === "log") { const b = S.grid.bounds_um; S.case.size_bins = { mode: "log", n: b.length - 1, dmin_um: b[0], dmax_um: b[b.length - 1] }; }
    else S.case.size_bins = { mode: "explicit", bounds_um: S.grid.bounds_um.map(x => +x.toPrecision(5)) };
    changed("size_bins"); render();
  } }, l)));
  const card = h("div", { class: "card" }, h("h2", {}, t("size_bins")), h("div", { class: "row" }, h("span", { class: "muted small" }, t("bins_mode")), seg));
  if (mode === "log") card.append(h("div", { class: "fields", style: { marginTop: "12px" } },
    field(t("n_bins"), num("size_bins.n", { step: 1 })), field(t("dmin"), num("size_bins.dmin_um")), field(t("dmax"), num("size_bins.dmax_um"))));
  else {
    const b = sb.bounds_um;
    card.append(h("div", { style: { marginTop: "12px" } }, editTable({ rows: b.map((_, i) => String(i)), cols: [t("bounds_um")], corner: "#",
      get: i => b[i], set: (i, _, v) => { b[i] = v; changed("size_bins"); }, onRemoveRow: i => { if (b.length > 2) { b.splice(i, 1); changed("size_bins"); render(); } } }),
    h("button", { class: "btn sm", style: { marginTop: "8px" }, onclick: () => { b.push(+(b[b.length - 1] * 2).toPrecision(4)); changed("size_bins"); render(); } }, "+ " + t("add_row"))));
  }
  // preview on a log axis
  const bd = S.grid.bounds_um, W = 760, H = 70, lo = Math.log10(bd[0]), hi = Math.log10(bd[bd.length - 1]);
  const X = v => 20 + (Math.log10(v) - lo) / (hi - lo || 1) * (W - 40);
  const svg = `<svg class="binsvg" viewBox="0 0 ${W} ${H}" width="100%">` + bd.slice(0, -1).map((v, i) =>
    `<rect x="${X(v)}" y="14" width="${Math.max(X(bd[i + 1]) - X(v) - 2, 1)}" height="26" rx="3"></rect><text x="${(X(v) + X(bd[i + 1])) / 2}" y="31" text-anchor="middle" font-size="10" fill="currentColor">${i + 1}</text>`).join("") +
    bd.map(v => `<text x="${X(v)}" y="58" text-anchor="middle" font-size="10" fill="currentColor" opacity=".7">${fmt(v)}</text>`).join("") + "</svg>";
  card.append(h("h3", {}, `${bd.length - 1} ${t("bins_n")} · µm`), h("div", { html: svg }));
  body.append(card);
}
