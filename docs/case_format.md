# SCRAM box app 运行文件格式与实现说明（v1.0，2026-10-08）

一个 JSON（case）描述一次 box 模拟。网页界面读写它，命令行 `python -m scrambox run case.json -o RUN`
也直接用它，用户不再接触 cfg。本文档描述已实现的格式；v0.2 草案中的设计均已落实，差异在各节注明。

## 1. 架构

```
网页界面 / CLI ──▶ case.json ──scrambox.compile──▶ RUN/ ─────────────▶ kernel/build/<v1|v23>/scram_box
                                    ├ case.resolved.json  默认值补全后的 case            │
                                    ├ grid.json           粒径段、类群、物种表           ▼
                                    ├ scram_species.dat   内核物种表                 output.txt ─▶ results.npz
                                    ├ run.nml             数值与开关                              summary.json
                                    └ inputs.txt          数组（关键字块）                        post.json（后处理）
```

- Fortran 驱动 `kernel/driver/scram_box.F90` 不解析 JSON，按内核的标准接口和顺序调用内核：
  每步 排放与气体源 → 稀释 → SOA 光解损失 exp(−0.0004·J_NO2·dt) → `scram_kernel_set_met` →
  `scram_kernel_set_hygro_rh(min(RH,0.98))` → `set_state` → n 次 `scram_kernel_step(dt/n)` → `get_state`。
- 内核来自 SCRAM 内核快照（ISORROPIA II v2.3 +
  短调用修正 07f62ee），`kernel/build_kernel.sh` 编出 v1 与 v23 两个可执行文件，`KERNEL_VERSION` 记录提交号。
  `kernel.isorropia = none` 用 v1 可执行文件并关闭 ISORROPIA（no-ISO 路径）。

## 2. 顶层字段

| 字段 | 界面位置 | 说明 |
|---|---|---|
| `case` | 情景 | name、description、scenario_id、category、reference、derived_from |
| `kernel` | 配置 › 内核与过程 | `isorropia`: `v23` / `v1` / `none` |
| `environment` | 配置 › 环境与时间 | T、P、RH、J_NO2、稀释、时间序列 |
| `time` | 配置 › 环境与时间 | duration_s、dt_s、**substep_s**、output_every_s |
| `species` | 配置 › 物种与组 | `custom`: 自定义惰性物种 |
| `mixing` | 配置 › 颗粒类群 | NC1 / NC6 / NC11 / custom |
| `size_bins` | 配置 › 粒径段 | 对数自动或手动边界 |
| `initial`、`emissions` | 初始与排放 | 颗粒（三种填法）、气体、时间廓线 |
| `processes`、`advanced` | 配置 › 内核与过程 | 过程开关与数值方案，默认 默认值 |
| `postprocess` | 后处理 | optics、ccn、hetchem |

## 3. 字段细节

### 3.1 时间与内核子步长

```json
"time": {"duration_s": 86400, "dt_s": 600, "substep_s": null, "output_every_s": 1800}
```

`substep_s` 用户可改；为 `null` 时取默认：**ISORROPIA 动态方案（表面蒸汽压，非平衡）10 s，
no-ISO 与整体平衡方案 600 s**（SCRAM `scram_target_dt`）。大于 dt_s 时截为 dt_s。

> 内核修正（2026-10-08）：`initstep` 先按剩余时间截断、再抬到 DTAEROMIN（60 s），
> 导致短于 60 s 的调用实际积分 60 s（纯凝并 1 h 数浓度损失：600/300/60 s 调用 1.77 %，10 s 为 9.4 %，
> 1 s 为 45 %）。现改为先下限后截断，
> 剩余时间 ≥ 60 s 的调用（全部 600 s 生产调用）结果逐位不变；修正后各子步长均为 1.77 %。

### 3.2 环境

```json
"environment": {
  "temperature_K": 298, "pressure_Pa": 101325, "relative_humidity": 0.7, "jno2_s1": 0,
  "time_series": {"time_s": [...], "temperature_K": [...], "relative_humidity": [...],
                  "pressure_Pa": [...], "jno2_s1": [...], "dilution_s1": [...]},
  "dilution": {"rate_s1": 1.5e-5, "background": "initial"}
}
```

时间序列线性插值、范围外取端点。稀释 x ← x − (1−e^(−k·dt))(x − x_bg)，背景默认为初始状态，
也可给 `{"particles": {...}, "gas_ugm3": {...}}`。

### 3.3 物种

物种表来自内核：ISO 版 14 个气溶胶物种（PMD PBC PNA PSO4 PNH4 PNO3 PHCL PSOA1P–5P PPOAlP PH2O，
有机为 CESM VBS），no-ISO 版 11 个；组顺序 BC DU SLT SUF POA SED（第 1 组为主组）。
自定义物种插在 PH2O 之前，为惰性物种：

```json
"species": {"custom": [{"name": "X1", "group": "SED", "density_gcm3": 1.5, "mw_gmol": 150, "kappa": 0.1}]}
```

### 3.4 类群与粒径段

类群生成规则同 `carma_model_mod`；NC1 = 内混，NC6 = 界限 0/0.5/1，NC11 = 0/0.5/0.9/1。
`custom` 时给 `external`、`primary_fraction_bounds`、`secondary_fraction_bounds`。
粒径段：`{"mode": "log", "n", "dmin_um", "dmax_um"}` 或 `{"mode": "explicit", "bounds_um": [...]}`。

### 3.5 初始颗粒与排放（三种填法）

| `mode` | 内容 | 类群归属 |
|---|---|---|
| `lognormal` | `modes`: name、number_cm3（排放为 number_cm3_s）、dg_um（数中值直径）、sigma_g、mass_fraction | 每个模态按干组成落到对应类群 |
| `table` | number_cm3[类群][粒径]，mass_fraction[类群]{物种} | 直接指定 |
| `bin_mass` | mass_ugm3{物种: [各粒径]}（排放为 mass_ugm3_s），number 可选（空则按中值粒径推算） | `mixing`: `internal` / `by_group` / `by_group_set`（`sets` 列表） |

气体：`initial.gas_ugm3`、`emissions.gas_ugm3_s`。时间廓线：
`emissions.time_profile`（颗粒排放倍数；为 null 时恒定）；`emissions.gas_time_profile`：
null = 同颗粒排放，`"jno2"` = 按 J_NO2/max(J_NO2)（光化学生成的替代），或 `{time_s, factor}`。
box 没有气相化学，H2SO4、HNO3、SOA 前体用气体源项给出。

### 3.6 在线 N2O5 非均相吸收（`chemistry.n2o5_uptake`）

```json
"chemistry": {"n2o5_uptake": {"enabled": true, "mode": "production",
                              "n2o5_production_ugm3_s": 1.2e-3, "k_other_s1": 2.0e-3,
                              "rule": "cesm", "gamma_cesm": 0.1, "species_gamma": null,
                              "clno2_yield": "bertram_thornton"}}
```

在 Fortran 驱动里每步（SCRAM 调用之前）计算：各粒径段 k_p = sfc_p/(0.5·dm_p/Dg + 4/(c·γ_p))，
表面与湿直径同 `scram_hetchem_sad_export`；`rule: cesm` 为 CESM 固定 γ 作用于 f_chm 表面，`species` 为物种 γ
（缺省取后处理 γ 表的 N2O5 列）。`mode: production` 时 [N2O5] = P/(Σk_p + k_other)（夜间 NO3+NO2 生成受限的稳态），
`concentration` 时为给定浓度（不消耗；只适合短时间敏感性试验，长时间会使硝酸盐–含水–表面正反馈失控）。
每吸收 1 mol N2O5：(2−φ) mol HNO3 进入气相，由内核热力学分配；φ mol ClNO2 释放并扣除颗粒氯，
φ = 1/(1 + [H2O]/(483[Cl⁻]))（Bertram & Thornton 2009，按各段颗粒含水与氯计算）。
输出 `het` 行：k、湿表面积、φ、累计 N2O5 吸收、HNO3 与 ClNO2 生成、[N2O5]。
在线 k 与后处理 cesm 规则的 k 逐位一致（同一状态下 4.98e-3 s⁻¹）。

## 4. 情景库

`scenarios/*.json` 由 `tools/make_scenarios.py` 生成（转换过程可复现）；用户情景存于
`SCRAMBOX_USER_DIR/scenarios`（默认 `~/.scrambox`）。界面为带搜索的分组下拉框，选中后显示信息面板
（环境、内核与过程、背景颗粒 N/PM1/PM2.5/PM10 与物种、背景气体、颗粒与气体排放、说明与文献），
面板数据由 `/api/describe` 从 case 本身计算，保证与实际输入一致。

| id | 类别 | 内容 |
|---|---|---|
| tutorial_minimal | 教学 | 新鲜 BC + 硫酸盐凝并 1 h，no-ISO |
| gmd_hazy_condensation / gmd_hazy_coag_cond | GMD 验证 | Seigneur hazy 三模态，一半硫酸盐一半 BC，H2SO4 2.29e-4 µg m⁻³ s⁻¹，12 h |
| gmd_paris_emission_only / coagulation / condensation / full | GMD 巴黎 A–D | 旧 BoxApp baseline12h：7 段初始质量与数浓度，BC/沙尘排放，气体源 |
| urban_plume_aging | 混合态 | PartMC urban plume：背景 Aitken/积聚模态，柴油/汽油/烹饪排放 12 h，稀释 1.5e-5 s⁻¹，按 J_NO2 的 H2SO4/HNO3/SOA 源，24 h |
| seasalt_nitrate | 混合态 | 海盐分配（PNA 0.3866、PHCL 0.6134）+ 大陆硫酸铵 + HNO3 |
| dust_crustal | 混合态 | 沙尘细/粗模态（CMAQ 地壳离子由 v2.3 内核诊断）+ 污染气体，RH 30→85 % |
| nocturnal_n2o5_clno2 | 混合态 | 沿海夜间 18–06 时：海盐、硫酸铵-硝酸铵、有机包覆颗粒，在线 N2O5 吸收（1 ppb/h），HNO3 与 ClNO2 |

旧情景转换要点（写在各情景 notes 中）：H2O 有机物按 298 K 的 C* 归入最近的 VBS 档
（C* > 1000 µg m⁻³ 的气相不表示）；颗粒 POA → PPOAlP（CESM POA 不挥发），POA 气体与排放按 C* 归入
SOA 气体档；旧的外混初始按 沙尘 / BC / 无机 / 有机 四个群分配，数浓度按体积分配。
海盐不采用 CMAQ sea_spc（含 SO4 而缺少 6 % 的 Mg/Ca/K，v1 下阴离子过剩，1 h 内约 25 % 的 Cl 以 HCl 逸出）。

> v2.3 海盐电荷修正：海盐以电中性 NaCl 进入 SCRAM（PNA:PHCL = 23:36.5），
> 原 v2.3 适配层把按 Na 诊断的 Ca/K/Mg 叠加上去，阳离子当量多出约 29 %（box 检验：各海盐段阳/阴离子比 1.22–1.25）。
> 现改为把钠当量按海水比例重新分配（Na 0.775、Mg 0.175、Ca 0.034、K 0.017），比值回到 0.98（余量为 Cl 摩尔质量记法差）。
> 只改 `USE_ISORROPIA2` 内的代码，默认预处理源码逐字节不变。12 h 后 HCl 释放：v1 3.2、原 v2.3 2.5、修正后 2.8 µg m⁻³。

## 5. 后处理（`scrambox/post`）

每个输出时刻、每个（粒径段，类群）为一个颗粒群；结果同时给出类群分辨与"每个粒径段合并为内混"的参考。

### 5.1 光学（`optics.py`，对应 `CARMA_DiagnoseBulk`）

```json
"postprocess": {"optics": {"enabled": true, "mixing": 2, "perclass": true,
                           "prognostic_moments": false, "dryshell_selector": 0, "every": 1}}
```

- mixing 0 按类型 Mie 表、1 体积平均折射率（在线 Mie）、2 简化壳 LUT（生产）、3 核壳 LUT；
  折射率表、物种→光学类型、κ 表、MAM4 κ-Köhler 湿半径（含 0.35/0.80 滞后）与内核逐项一致；
  LUT 用 `kernel/fetch_luts.sh` 从内核运行目录复制（波长与 CARMA/RRTMG 网格一致，SW 第 10 波段 0.533 µm）。
- 输出：14 个短波波段消光、SSA、g，16 个长波吸收，550 nm 消光/SSA/g、BC 与沙尘贡献，按类群与粒径段的消光，
  BC 吸收增强 E_abs（相对同体积裸 BC 球）。box 中 τ 换成消光系数（Mm⁻¹）。
- 刻意保留的行为（非生产路径）：mixing 0 的 RH 节点按均匀网格索引而 mie_rh 非均匀
  （RH 0.82 读到 0.95–0.98 的表，吸湿增长高估）；mixing 1 的混合折射率不含水却用于湿半径。

### 5.2 云活化（`ccn.py`，对应 `scram_cloud_aerosol_export`）

```json
"postprocess": {"ccn": {"enabled": true, "supersaturations_pct": [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0],
                        "updrafts_m_s": [0.2, 0.5, 1.0, 2.0]}}
```

每个群取粒径段参考干半径与体积加权 κ，sc = 2A·√(A/(27κr³))；粒径段内按由 N 与干体积反推的线性质量倾斜
计算部分活化；给定 S 的 CCN 谱（CESM CCN3 路径），以及 Abdul-Razzak & Ghan 最大过饱和度（f1=0.5、f2=1，
ndrop_bam 热力学，w ≥ 0.2 m/s）下的云滴数、各类群活化比例。

### 5.3 非均相化学（`hetchem.py`）

k = Σ_p sfc_p / (0.5·dm_p/Dg + 4/(c·γ_p))（CAM `hetrxtrate`，Dg = 0.1 cm² s⁻¹，c = c0·√T）。
表面积与湿直径同 SCRAM `scram_hetchem_sad_export`（干半径由 N 与 M 求得，湿直径用本群 PH2O）。

| 规则 | γ |
|---|---|
| `cesm` | CESM 固定值（N2O5 0.1、HO2 0.2、NO3 1e-3、NO2 1e-4），作用于 f_chm 部分的表面 |
| `species` | 用户可编辑的物种 × 气体 γ 表，按干体积分数加权，作用于全部湿表面；默认表复现 CESM |
| `bertram_thornton` | N2O5：Bertram & Thornton (2009)，由本群含水、硝酸根、氯离子浓度计算；可选有机包覆限制（Anttila 2006 / Riemer 2009） |

输出 k、有效 γ（以总湿表面计）、寿命，类群分辨与内混参考对照。

## 6. 实验组（`experiments.py`）

```json
{"name": "...", "base": {case}, "postprocess": {...},
 "variants": [{"label": "NC1", "overrides": {"mixing.preset": "NC1"}}, ...]}
```

变体并行运行（`SCRAMBOX_EXP_WORKERS`，默认 4），对比指标：N、干 PM、PM1、含水、各物种、550 nm 消光
（类群分辨与内混）、SSA、E_abs、各 S 的 CCN、各 w 的 Nd、各气体各规则的 k_het。

## 7. 运行方式

```bash
python3 -m scrambox serve --port 8765          # 服务器上每个用户各起一个实例
ssh -L 8765:localhost:8765 <server>              # 本地浏览器打开 http://localhost:8765
python3 -m scrambox run scenarios/dust_crustal.json -o runs/dust   # 命令行
```

## 8. 状态与后续

已完成：内核驱动、编译器、11 个情景、三个后处理模块、在线 N2O5 吸收（HNO3 / ClNO2 回到热力学）、
实验组、网页界面（中英文）。
后续：单机安装包（内核静态链接 + 打包 Python 与网页）
（需要从全球模式输出导出单柱状态）；PartMC 对比（同一情景输入）。
