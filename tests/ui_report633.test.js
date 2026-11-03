// Behavioral regressions for report 633. Run with node ui_report633.test.js <web-dir>.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');

const dir = process.argv[2] || __dirname;
class Element {
  constructor(tag, text = '') { this.tag = tag; this.children = []; this.attrs = {}; this.textContent = text; this.nodeType = 1; this.style = {}; this.classList = { add() {}, remove() {} }; }
  append(...children) { this.children.push(...children); }
  setAttribute(k, v) { this.attrs[k] = v; if (k === 'value') this.value = v; }
  addEventListener(k, fn) { this['on' + k] = fn; }
  click() { if (this.onclick) return this.onclick({ target: this }); }
  get innerHTML() { return this._html || ''; }
  set innerHTML(v) { this._html = v; this.children = []; }
}
const chip = new Element('span');
const context = vm.createContext({
  console, setTimeout, clearTimeout, setInterval: () => 0, Charts: { fmt: String },
  URL: { createObjectURL: () => 'blob:test', revokeObjectURL() {} }, URLSearchParams, Blob,
  location: { search: '' }, localStorage: { getItem: () => null },
  document: {
    createElement: tag => new Element(tag), createTextNode: text => ({ nodeType: 3, textContent: String(text) }),
    querySelector: sel => sel === '#kernelChip' ? chip : new Element('div'),
    querySelectorAll: () => [], body: new Element('body')
  },
});
function load(name) {
  let src = fs.readFileSync(path.join(dir, name), 'utf8');
  if (name === 'app_results.js') src = src.replace(/init\(\)\.catch\([\s\S]*?\);\s*$/, '');
  vm.runInContext(src, context, { filename: name });
}
load('i18n.js'); load('app.js'); load('app_inputs.js'); load('app_results.js');
vm.runInContext('S.lang = "en"; render = () => {}; changed = () => {};', context);
const S = vm.runInContext('S', context);
const evaluate = source => vm.runInContext(source, context);
context.originalApi = evaluate('api');
const find = (root, tag, text) => {
  const all = [];
  function walk(n) { if (n && n.tag) { all.push(n); n.children.forEach(walk); } }
  walk(root);
  return all.find(n => n.tag === tag && (text === undefined || n.children.some(c => c.textContent === text)));
};
const baseCase = () => ({
  case: { name: 'GMD D' }, kernel: { isorropia: 'v23' }, mixing: { preset: 'NC6' },
  size_bins: { mode: 'explicit', bounds_um: [0.00316227766016838, 0.01, 0.1] },
  initial: { particles: { mode: 'bin_mass', mixing: 'by_group_set', sets: [['PBC'], ['PSO4']],
    number_cm3: [12, 13], mass_ugm3: { PBC: [1, 2], PSO4: [3, 4] } } },
  emissions: { particles: { mode: 'lognormal', modes: [] } },
});
S.case = baseCase();
S.grid = { bounds_um: [0.00316227766016838, 0.01, 0.1], class_labels: ['A', 'B'], aerosol: ['PBC', 'PSO4', 'PH2O'] };

const failures = [];
async function test(name, fn) { try { await fn(); process.stdout.write(`PASS ${name}\n`); } catch (e) { failures.push(name); process.stdout.write(`FAIL ${name}: ${e.message}\n`); } }
(async () => {
  await test('B1 input mode round trip retains mass, number and mixing sets', () => {
    let card = evaluate('particleEditor("initial", false)');
    find(card, 'button', 'Log-normal modes').click();
    card = evaluate('particleEditor("initial", false)');
    find(card, 'button', 'Mass per species and bin').click();
    assert.equal(JSON.stringify(S.case.initial.particles), JSON.stringify(baseCase().initial.particles));
  });
  await test('B3 manual bounds survive Auto and Manual round trip exactly', () => {
    S.case.size_bins = baseCase().size_bins;
    S.manualBounds = null;
    let body = new Element('div'); evaluate('setupBins')(body);
    find(body, 'button', 'Log-spaced (auto)').click();
    body = new Element('div'); evaluate('setupBins')(body);
    find(body, 'button', 'Manual bounds').click();
    assert.equal(JSON.stringify(S.case.size_bins.bounds_um), JSON.stringify(baseCase().size_bins.bounds_um));
  });
  await test('B4 by_group_set is the selected visible mixing option', () => {
    S.case.initial.particles = baseCase().initial.particles;
    const card = evaluate('particleEditor("initial", false)');
    const selected = find(card, 'option', 'by listed species sets');
    assert.ok(selected && Object.hasOwn(selected.attrs, 'selected'));
    assert.ok(JSON.stringify(S.case.initial.particles).includes('by_group_set'));
  });
  await test('B4 submitted mixing and species sets are echoed in the run summary', () => {
    S.model = 'scram'; S.meta = { partmc: {}, limits: {} }; S.runs = []; S.runId = null;
    const page = new Element('div'); evaluate('pageRun')(page);
    const tree = JSON.stringify(page);
    assert.ok(tree.includes('by_group_set'));
    assert.ok(tree.includes('1=[PBC]'));
    assert.ok(tree.includes('2=[PSO4]'));
  });
  await test('B2 unit guidance states the cm^-3 to m^-3 conversion without changing values', () => {
    const before = JSON.stringify(S.case.initial.particles.number_cm3);
    const card = evaluate('particleEditor("initial", false)');
    assert.ok(JSON.stringify(card).includes('10⁶ m⁻³'));
    assert.equal(JSON.stringify(S.case.initial.particles.number_cm3), before);
  });
  await test('B6 PartMC run is disabled for bin_mass input with an explanation', () => {
    S.model = 'partmc'; S.meta = { partmc: { available: { v23: true }, defaults: {} }, limits: {} };
    S.runs = []; S.runId = null;
    const page = new Element('div'); evaluate('pageRun')(page);
    const tree = JSON.stringify(page);
    assert.ok(tree.includes('PartMC currently accepts log-normal'));
    assert.ok(tree.includes('"disabled":""'));
    S.model = 'scram';
  });
  await test('B8 selected scenario is visibly distinguished from editor case', () => {
    S.meta = { categories: { study: { en: 'Study', zh: '研究' } } };
    S.scenario = { id: 'selected', name: { en: 'Selected' }, summary: { en: 'Summary' }, category: 'study',
      builtin: true, tags: [], case: { case: {} } };
    S.scenarios = [S.scenario]; S.describe = { error: 'skip' }; S.loadedScenarioId = 'other';
    const page = new Element('div'); evaluate('pageScenario')(page);
    assert.ok(JSON.stringify(page).includes('Selected, not loaded'));
  });
  await test('B9 advanced controls have visible names and meaningful unit hints', () => {
    S.case = baseCase();
    S.case.advanced = { bc_core_floor_um: 0.05, salt_core_floor_um: 0.05, dust_core_floor_um: 0.9,
      nucl_rate_cap_cm3_s: 1000, nucl_h2so4_gas_fraction_max: 0.25, nucl_nh3_gas_fraction_max: 0.25,
      salt_rdb_bidirectional: false, salt_fix_size_rdb: false, dust_fix_size_rdb: false };
    S.case.processes = { redistribution: 'moving_center_dualpivot', nucleation: { enabled: false }, density: { fixed_kgm3: 1400 } };
    S.case.chemistry = { n2o5_uptake: { enabled: false } };
    S.meta = { kernels: { v23: { available: true, version: '' }, v1: { available: true, version: '' } },
      enums: { 'processes.redistribution': ['moving_center_dualpivot'] }, defaults: { advanced: S.case.advanced,
        processes: { density: { fixed_kgm3: 1400 } } } };
    const body = new Element('div'); evaluate('setupKernel')(body);
    const tree = JSON.stringify(body);
    for (const label of ['BC core minimum diameter (µm)', 'Nucleation rate cap (cm⁻³ s⁻¹)',
      'Fixed particle density (kg m⁻³)', 'Salt redistribution in both directions']) assert.ok(tree.includes(label), label);
  });
  await test('B5 invalid results response shows an error instead of rejecting', async () => {
    S.runId = 'run-1'; S.runs = [{ id: 'run-1', state: 'done' }]; S.results = null;
    evaluate('api = async () => { throw new Error("invalid JSON") }');
    const box = new Element('div');
    try { await evaluate('showRun')(box); }
    finally { evaluate('api = originalApi'); }
    assert.ok(JSON.stringify(box).includes('invalid JSON'));
  });
  await test('B5 malformed HTTP 200 JSON is treated as an API error', async () => {
    context.fetch = async () => ({ status: 200, ok: true, statusText: 'OK', json: async () => { throw new SyntaxError('bad JSON'); } });
    await assert.rejects(evaluate('api("runs/run-1/results")'), /invalid JSON/);
  });
  await test('B11 open-system mass change is visibly inapplicable without a pass threshold', () => {
    // Synthetic B5 run summary plus conservation_diagnostic from its resolved case.
    const summary = JSON.parse(fs.readFileSync(path.join(__dirname, 'fixtures', 'report633_conservation_summary.json'), 'utf8'));
    const card = evaluate('conservationCard')(summary);
    const tree = JSON.stringify(card);
    assert.ok(tree.includes('Gas + particle mass change'));
    assert.ok(tree.includes('Not applicable: particle emissions, gas sources'));
    assert.ok(tree.includes('PSOA3P+SOA3P'));
    assert.ok(tree.includes('22.915'));
    assert.ok(tree.includes('PSOA5P+SOA5P'));
    assert.ok(tree.includes('7.268'));
    assert.ok(tree.includes('dimensionless'));
    assert.ok(tree.includes('no pass threshold is defined'));
    assert.ok(!tree.includes('PASS'));
    assert.ok(!tree.includes('NaN'));
  });
  await test('B11 closed-box diagnostic is identified without declaring success', () => {
    const card = evaluate('conservationCard')({ conservation_rel: { 'PSO4+SULF': 1 / 3 },
      conservation_diagnostic: { definition: 'Initial-to-final gas and particle mass relative change',
        units: 'dimensionless', applicable: true, reasons: [], pass_threshold: null } });
    const tree = JSON.stringify(card);
    assert.ok(tree.includes('Closed-box diagnostic applies'));
    assert.ok(tree.includes('PSO4+SULF'));
    assert.ok(!tree.includes('PASS'));
  });
  await test('B7 PartMC progress shows observed files and timing in list and detail', () => {
    const run = { id: 'pm-1', name: 'Synthetic', model: 'partmc', state: 'running', progress: 0.04,
      output_count: 3, expected_output_count: 72, elapsed_s: 120,
      last_output_at: 1000, last_output_age_s: 36, log_bytes: 1234 };
    S.runs = [run]; S.runId = run.id;
    const list = evaluate('runList')();
    assert.ok(JSON.stringify(list).includes('3 / 72 output files'));
    const box = new Element('div');
    evaluate('showRun')(box);
    const detail = JSON.stringify(box);
    assert.ok(detail.includes('last output 36 s ago'));
    assert.ok(detail.includes('log 1234 B'));
    assert.ok(!detail.includes('stalled'));
  });
  if (failures.length) process.exitCode = 1;
})();
