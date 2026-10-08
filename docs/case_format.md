# SCRAM box app 运行文件格式（草案 v0.1，2026-10-08）

一个 JSON 文件描述一次 box 模拟。GUI 读写它，命令行也能直接用它，用户不再接触 cfg。
字段名用英文，界面上显示中文。本文是给讨论用的初稿，确认后再写代码。

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

## 4. 需要确认的问题

1. 用 Python 编译器展开、Fortran 只读展开后数组的架构可以吗？
2. 初始颗粒默认用"对数正态模态"填法、表格作为高级选项，是否合适？
3. 气体源项（不含气相化学）作为 box 的标准做法可以吗？以后如果要接 CAMP 气相化学再扩展。
4. 现有 7 个典型场景（教学、hazy、巴黎 A–D）的原始数值在旧 BoxApp 的模板里，我会逐个转写成新格式。
