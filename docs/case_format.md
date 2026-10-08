# SCRAM box app 运行文件格式（草案 v0.2，2026-10-08）

一个 JSON 文件描述一次 box 模拟。GUI 读写它，命令行也能直接用它，用户不再接触 cfg。
字段名用英文，界面上显示中文。

v0.1 的四点已确认（Python 编译器展开、模态为默认初始填法、气体源项代替气相化学、转写 7 个旧场景）。
v0.2 新增：情景库与下拉选择（第 5 节）、后处理模块——光学辐射、云激活、非均相化学吸收（第 6 节）、
多配置对比实验（第 7 节）、开发顺序（第 8 节）。

## 1. 总体架构

```
GUI (浏览器)  ──编辑──▶  case.json  ──Python 编译器──▶  run/ 目录  ──▶  scram_box (Fortran)  ──▶  结果
                                       ├─ run.nml            数值和开关
                                       ├─ scram_species.dat  物种表（内核现有解析器读取）
                                       └─ inputs.nc          初始状态、排放、气象时间序列
```

- **JSON 是唯一的用户格式**；`run.nml`、`scram_species.dat` 是生成的中间文件，用户不需要看。
- **Fortran 端不解析 JSON**：由 Python 校验并展开（类群生成、粒径分档、分布映射到格点），Fortran 只读展开后的数组。这样内核可执行文件不依赖 JSON 库，单机安装包更容易做。
- 内核走与 CESM **完全相同的接口**：`scram_cesm_init_from_discretization` → `scram_cesm_set_met` →
  `scram_cesm_set_state` → 循环 `scram_cesm_step` → `scram_cesm_get_state`。
- ISORROPIA 版本对应两个内核可执行文件（v1、v2.3），由 `kernel.isorropia` 选择；`none` 用 v1 文件，
  `use_isorropia = 0` 走 no-ISO 路径。

## 2. 顶层结构

| 字段 | 对应 GUI 标签页 | 说明 |
|---|---|---|
| `schema` | — | `"scram-box/0.1"` |
| `case` | 情景 | 名称、描述、来源模板 |
| `kernel` | 顶栏 | ISORROPIA 版本、阳离子分配谱 |
| `environment` | 情景 | 温度、气压、RH（常数或时间序列） |
| `time` | 情景 | 时长、步长、输出间隔 |
| `species` | 化学组分 | 气溶胶物种、气体物种，每个气溶胶物种归属一个组 |
| `mixing` | 颗粒类群 | 内混或外混；外混时的主组、分段边界；生成的类群列表（只读回显） |
| `size_bins` | 粒径分档 | 对数自动或手动边界 |
| `initial` | 初始与排放 | 初始颗粒（模态或表格）和初始气体 |
| `emissions` | 初始与排放 | 颗粒排放、气体源（含 H2SO4、SOAG0–4 生成率） |
| `processes` | 过程与数值 | 凝并、冷凝、成核、热力学方案、求解器、再分配 |
| `advanced` | 过程与数值的"高级"折叠区 | 核心粒径下限、成核限幅等，默认取 CESM 生产值 |
| `output` | 运行与结果 | 输出变量和格式 |

## 3. 各部分

### 3.1 `kernel`

```json
"kernel": {
  "isorropia": "v23",                       // "v23" | "v1" | "none"
  "crustal_speciation": "cmaq_aero6_7"      // 仅 v23 有效；也可 "user" 并给出 crustal_table
}
```

### 3.2 `environment` 和 `time`

```json
"environment": {
  "temperature_K": 298.0,
  "pressure_Pa": 101325.0,
  "relative_humidity": 0.70,
  "time_series": null          // 或 {"time_s": [...], "temperature_K": [...], "relative_humidity": [...]}
},
"time": { "duration_s": 43200, "dt_s": 600, "output_every_s": 600 }
```

### 3.3 `species`

默认从物种库载入当前 hetfix 的 14 个气溶胶物种（含 PNH4/PNO3/PHCL，有机为 VBS 五档）和全部气体。
用户可以加自定义物种，只需选组、填密度和分子量；其余热力学参数在"高级属性"里，有默认值。

```json
"species": {
  "library": "scram_hetfix_iso14",
  "aerosol": [
    {"name": "PMD",    "group": "DU",  "density_gcm3": 2.60, "mw_gmol": 28.0},
    {"name": "PBC",    "group": "BC",  "density_gcm3": 1.70, "mw_gmol": 12.0},
    {"name": "PNA",    "group": "SLT", "density_gcm3": 1.90, "mw_gmol": 23.0},
    {"name": "PSO4",   "group": "SUF", "density_gcm3": 1.77, "mw_gmol": 98.0},
    {"name": "PSOA1P", "group": "SED", "density_gcm3": 1.00, "mw_gmol": 250.0}
  ],
  "groups": ["BC", "DU", "SLT", "SUF", "POA", "SED"],
  "gas": "library"              // 气体列表按物种库；用户一般不改
}
```

### 3.4 `mixing`（颗粒类群）

```json
"mixing": {
  "external": true,
  "primary_fraction_bounds":   [0.0, 0.5, 0.9, 1.0],
  "secondary_fraction_bounds": [0.0, 0.5, 0.9, 1.0],
  "preset": "NC11"              // "NC1" | "NC6" | "NC11" | "custom"
}
```

生成规则与 CESM 完全相同（`carma_model_mod`）：第 1 组为主组，中间各组为次组，最后一组为余量；
保留各组下限之和小于 1 的组合。GUI 实时显示生成的类群数和每个类群的组成区间。

### 3.5 `size_bins`

```json
"size_bins": { "mode": "log", "n": 7, "dmin_um": 0.003, "dmax_um": 10.0 }
// 或
"size_bins": { "mode": "explicit", "bounds_um": [0.003, 0.01, 0.0398, 0.1585, 0.631, 2.5, 5.0, 10.0] }
```

### 3.6 `initial`

两种填法，GUI 用下拉切换：

```json
"initial": {
  "particles": {
    "mode": "lognormal",
    "modes": [
      {"name": "Aitken",       "number_cm3": 8000, "dg_um": 0.03, "sigma_g": 1.6,
       "mass_fraction": {"PSO4": 0.7, "PSOA1P": 0.3}},
      {"name": "accumulation", "number_cm3": 1500, "dg_um": 0.15, "sigma_g": 1.5,
       "mass_fraction": {"PSO4": 0.5, "PNH4": 0.15, "PBC": 0.1, "PSOA2P": 0.25}}
    ]
  },
  "gas_ugm3": {"SULF": 0.0, "NH3": 5.0, "HNO3": 2.0, "HCL": 0.1, "SOA1P": 0.5}
}
```

或者表格模式（行列名自动生成，只填数字）：

```json
"particles": {
  "mode": "table",
  "number_cm3":   [[...每个类群一行，每个粒径档一列...]],
  "mass_fraction": {"by": "class", "table": [[...每个类群一行，每个物种一列...]]}
}
```

模态模式下，Python 编译器把每个模态积分到粒径档，再按组成落到对应类群；表格模式直接使用。

### 3.7 `emissions`

```json
"emissions": {
  "particles": { "mode": "lognormal", "modes": [ {"number_cm3_s": 1.0, "dg_um": 0.05, "sigma_g": 1.8,
                                                   "mass_fraction": {"PBC": 0.4, "PPOAlP": 0.6}} ] },
  "gas_ugm3_s": {"SULF": 1.0e-4, "NH3": 2.0e-4, "SOA3P": 5.0e-5},
  "time_profile": null      // 或分时段倍数
}
```

气体源项让 box 不需要气相化学，也能与 PartMC 用同一组预设输入对比。

### 3.8 `processes`

```json
"processes": {
  "coagulation": true,
  "coagulation_repartition": "nearest",      // "nearest"（解析，无需系数文件）| "legacy"
  "condensation": true,
  "nucleation": {"enabled": false, "model": "binary"},
  "sulfate": "dynamic",                      // Sulfate_computation
  "thermodynamics": "dynamic",               // "dynamic"（表面蒸气压）| "equilibrium"（整体）
  "dynamic_solver": "etr",                   // "euler" | "etr" | "ros2"
  "redistribution": "moving_center_dualpivot",
  "wet_diameter": "isorropia",               // "isorropia" | "gerber"
  "density": {"mode": "variable", "fixed_kgm3": 1400}
}
```

默认值取当前 CESM 生产配置（NC6/NC11）。

### 3.9 `output`

```json
"output": { "format": "netcdf", "variables": ["number", "mass", "gas", "water", "ph", "wet_diameter"] }
```

结果写成一个 netCDF（时间 × 类群 × 粒径档 × 物种），外加一个运行摘要 JSON（耗时、守恒误差、
内核版本 `KERNEL_VERSION`、ISORROPIA 清单）。GUI 的结果页和对比功能都读这两个文件。

## 4. 已确认（2026-10-08）

1. Python 编译器展开 JSON，Fortran 只读展开后的数组。
2. 初始颗粒默认用对数正态模态填法，表格作为高级选项。
3. box 用气体源项，不做气相化学；以后需要再接 CAMP。
4. 转写旧 BoxApp 的 7 个典型场景。

## 5. 情景库

- **下拉菜单选择，不平铺卡片**。下拉框带搜索，按类别分组（教学、验证、城市、混合态、自定义）。
- 选中后，旁边的信息面板显示该情景的：环境条件、背景气体浓度、背景颗粒（总数浓度、质量、主要组分）、
  颗粒和气体排放（速率与时间变化）、文献出处。面板只读，点"载入"后才写入当前运行文件，各标签页可再改。
- **自定义情景**：任何当前运行文件都可以"另存为情景"，存入用户情景库（服务器版存在用户目录，
  单机版存在本地），之后出现在下拉菜单的"自定义"类别里，可重命名、删除、导出为 JSON 分享。
- 内置情景：教学（BC + 硫酸盐）、GMD hazy 冷凝、hazy 凝并 + 冷凝、巴黎 A–D（Zhu et al. 2015），
  以及计划书需要的城市羽流老化、海盐 + 硝酸盐、沙尘 + 地壳离子（v2.3）。

```json
"case": { "name": "城市羽流老化", "scenario_id": "mixing_urban_plume", "description": "...",
          "reference": "...", "derived_from": "builtin" }
```

### 5.1 计划书需要的三个混合态情景（用户确认必须有）

三者都对应计划书 WP2 的老化情景，并作为 PartMC 对比的共同输入（同一组背景、排放和气体源项）。
数值在实现时从下列来源取，并在信息面板里注明出处；若有组内观测（SP2、TEM-EDX），优先用观测约束。

| 情景 | 考察的过程 | 关键设定 | 参数来源（初定） |
|---|---|---|---|
| 城市羽流老化 `mixing_urban_plume` | 新鲜 BC/POA 被硫酸盐、硝酸铵、SOA 包覆，混合态从外混走向内混；N₂O₅、HO₂ 吸收随包覆变化 | 背景 + 交通/燃烧排放；H₂SO₄、HNO₃、NH₃、SOAG0–4 气体源项；日变化温湿度；24–48 h | PartMC urban plume 基准情景（Riemer et al., 2009 一系列）的排放与背景，按 NCP/YRD 观测调整 |
| 海盐 + 硝酸盐 `mixing_seasalt_nitrate` | HNO₃ 与海盐反应造成氯亏损，Cl⁻ 可得性决定 ClNO₂ 产率 φ | 沿海受污染气团：海盐粗粒子 + 城市细粒子；HNO₃、HCl、NH₃；高 RH | 海盐谱用 CMAQ `sea_spc`；气体与颗粒负荷参考沿海站点 N₂O₅/ClNO₂ 观测（计划书 WP1 的合作数据） |
| 沙尘 + 地壳离子 `mixing_dust_crustal` | 沙尘吸收 HNO₃/SO₂ 生成 Ca(NO₃)₂ 等吸湿盐，粗模态硝酸盐和含水量；**专门检验 v2.3 地壳离子与反向模式电荷平衡** | 沙尘粗细模态 + 污染细粒子；HNO₃、SO₂→H₂SO₄、NH₃；ISORROPIA v2.3 | 沙尘谱与地壳阳离子用 CMAQ `dust_spc`/`asoil_*_fac`（已入库）；负荷参考春季沙尘过境观测 |

三个情景都会给出 v1 与 v2.3 两套运行结果；沙尘情景在 v1 下 Ca/K/Mg 不起作用，可作为地壳离子效应的对照。

## 6. 后处理模块（与 CESM 同一套计算）

后处理读取模拟输出（每个输出时刻、每个类群 × 粒径档的数浓度、各物种质量、含水量、湿直径），
不改变模拟本身，因此可以对同一次模拟反复换参数分析。

### 6.1 光学与辐射

- **直接复用 CESM 的 SCRAM 光学模块**（`scram_reducedshell_optics_mod`、`scram_coreshell_optics_mod`，
  只依赖 netCDF，宿主模块写桩即可），以及同一份查找表（`scram_reducedshell_optics_lut_cdf5.nc`、
  `scram_coreshell_optics_lut.nc`）。
- 混合方式与 CESM 的 `Optics_mixing` 一致：0 按物种外混，1 体积平均折射率 Mie，2 简化壳层查找表，
  3 核壳查找表；`Optics_perclass`（按类群分别算或先合并）、`Optics_dryshell_selector` 也可选。
- 波段与 CESM 相同（RRTMG 短波波段）外加 550 nm 诊断；湿直径按 CESM 约定用环境 RH 和 κ 计算。
- 输出：消光、单次散射反照率、不对称因子（按类群、按粒径档、合计），以及吸收增强 E_abs。

```json
"postprocess": { "optics": { "enabled": true, "mixing": 2, "perclass": true,
                             "dryshell_selector": 2, "wavelengths": "rrtmg_sw+550nm" } }
```

### 6.2 云激活

- 与 CESM 的 `scram_cloud_aerosol_export` 同一算法：每个类群单独算 κ → Köhler 临界过饱和度 → 临界直径，
  粒径档内按单分散处理活化比例；可选 Abdul-Razzak & Ghan (2000) 上升速度求最大过饱和度。
- 该算法在 CESM 里嵌在 CAM 相关代码中，box 里按相同公式实现，并用 CESM 的诊断输出核对。
- 输出：给定过饱和度序列（默认 0.05、0.1、0.2、0.3、0.5、1.0 %）的 CCN 谱，按类群分解；
  给定上升速度时的最大过饱和度和活化数。

```json
"postprocess": { "ccn": { "enabled": true, "supersaturation_pct": [0.05, 0.1, 0.2, 0.3, 0.5, 1.0],
                          "updraft_m_s": [0.1, 0.5, 1.0] } }
```

### 6.3 非均相化学（吸收系数实验）

对应计划书 H1：同一颗粒群下，比较"按类群分别算 γ"与"按平均成分算 γ"与"CESM 固定 γ"。

- **γ 表**：行是物种（或类群），列是气体（N₂O₅、HO₂、NO₃、NO₂，可扩展），只填数字；
  另可选参数化（如 γ(N₂O₅) 的 Bertram & Thornton 2009、γ(HO₂) 的 Song 2020），参数化需要的
  含水量、Cu/Fe、有机包覆从模拟输出取。
- **混合规则**（下拉）：按类群分别计算再求和（类群分辨）、按体积/质量加权得到平均颗粒的 γ（内混）、
  CESM 固定值（γ(N₂O₅)=0.1、γ(HO₂)=0.2，作用于总表面积）。
- **计算**：每个类群 × 粒径档的湿表面积，一阶损失速率 k = Σ γ v A / 4（含 Fuchs–Sutugin 气相扩散修正），
  以及群体有效 γ = k / (v A_total / 4)。
- **输出**：k_het、群体有效 γ、类群分辨与内混之比（γ_pop / γ_mean，即计划书的判据量），随时间变化。
- 第一阶段只做诊断（不回馈颗粒成分）；第二阶段可选在线耦合：给定 N₂O₅ 生成率，吸收产生的 HNO₃
  （或按 φ 生成 ClNO₂）回到热力学，影响硝酸盐分配。

```json
"postprocess": { "hetchem": {
  "enabled": true,
  "gases": ["N2O5", "HO2"],
  "gamma_table": { "by": "species", "rows": ["PSO4", "PNH4", "PNO3", "PSOA1P", "PBC", "PMD", "PNA"],
                   "N2O5": [0.02, 0.02, 0.005, 0.001, 0.0, 0.01, 0.03],
                   "HO2":  [0.05, 0.05, 0.05, 0.01, 0.0, 0.05, 0.10] },
  "mixing_rules": ["class_resolved", "internal_mean", "cesm_fixed"],
  "online_coupling": false } }
```

（上例 γ 数值只是格式示意的占位数，不是推荐值；内置的默认 γ 表会从计划书的参数汇编
`parameter_compilation.xlsx` 取值并注明出处。）

## 7. 多配置对比实验

- 一个"实验组"由一个基准运行文件加若干变体组成，每个变体只写与基准不同的字段，例如：
  NC1 / NC6 / NC11、ISORROPIA v1 / v2.3、不同 γ 表、不同光学混合方式。
- 结果页并排比较：数浓度谱、各类群组成、含水量和 pH、光学量、CCN 谱、k_het 与 γ_pop/γ_mean。
- 后处理参数（光学方式、过饱和度、γ 表）的变体不需要重跑模拟，直接对已有输出重新分析。

```json
{ "schema": "scram-box-experiment/0.1", "base": "case.json",
  "variants": [ {"name": "NC1",  "set": {"mixing.preset": "NC1"}},
                {"name": "NC11", "set": {"mixing.preset": "NC11"}},
                {"name": "v1",   "set": {"kernel.isorropia": "v1"}} ] }
```

## 8. 开发顺序

1. 内核 JSON 驱动（Python 编译器 + Fortran runner），用旧 BoxApp 的教学和 hazy 场景做回归；
   接着测 ISORROPIA v2.3 和沙尘情景。
2. 后处理：光学（复用 CESM 模块）→ 云激活 → 非均相化学诊断。
3. 实验组与对比。
4. 网页界面：情景下拉与信息面板、各标签页、表格编辑器、结果与对比页。
5. 非均相化学在线耦合；单机安装包。
