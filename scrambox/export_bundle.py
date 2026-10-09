"""Build a self-contained, data-derived SCRAM Box result archive."""
from __future__ import annotations

import io
import json
import csv
import threading
import zipfile
from datetime import datetime, timezone

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


MAX_RESULTS_JSON = 96 * 1024 * 1024
_PLOT_LOCK = threading.Lock()


def _default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def _json_bytes(value):
    return json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False,
                      default=_default).encode("utf-8")


def _array(value, *, ndim=None):
    array = np.asarray(value, dtype=float)
    if ndim is not None and array.ndim != ndim:
        raise ValueError("result array has an unexpected shape")
    if not np.isfinite(array).all():
        raise ValueError("result array contains non-finite values")
    return array


def _picture(title, xlabel, ylabel, series, *, logx=False, logy=False,
             stacked=False, bar_labels=None):
    """Render actual numeric series; never emit a blank placeholder image."""
    if not series:
        raise ValueError(f"no numeric series for {title}")
    with _PLOT_LOCK:
        figure = Figure(figsize=(10, 5.5), dpi=140, layout="constrained")
        axis = figure.add_subplot(111)
        drawn = 0
        bottom = None
        for label, xs, ys in series:
            y = _array(ys, ndim=1)
            if not len(y):
                continue
            if bar_labels is not None:
                if len(y) != len(bar_labels):
                    raise ValueError(f"bar series length mismatch for {title}")
                if bottom is None:
                    bottom = np.zeros(len(y))
                axis.bar(np.arange(len(y)), y, bottom=bottom if stacked else None, label=label)
                if stacked:
                    bottom += y
            else:
                x = _array(xs, ndim=1)
                if len(x) != len(y):
                    raise ValueError(f"line series length mismatch for {title}")
                axis.plot(x, y, label=label, linewidth=1.8)
            drawn += 1
        if not drawn:
            raise ValueError(f"no numeric series for {title}")
        axis.set_title(title)
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.grid(True, alpha=0.25)
        if logx:
            axis.set_xscale("log")
        if logy:
            axis.set_yscale("symlog", linthresh=1e-10)
        if bar_labels is not None:
            axis.set_xticks(np.arange(len(bar_labels)), bar_labels)
        axis.legend(loc="best", fontsize=8)
        content = io.BytesIO()
        FigureCanvasAgg(figure).print_png(content)
        image = content.getvalue()
        if not image.startswith(b"\x89PNG\r\n\x1a\n") or len(image) < 2000:
            raise ValueError(f"invalid rendered image for {title}")
        return image


def _base_charts(result):
    time_h = _array(result["time_s"], ndim=1) / 3600.0
    labels = list(result["class_labels"])
    species = list(result["species"])
    number = _array(result["number_m3"], ndim=3) / 1e6
    mass = _array(result["mass_ugm3"], ndim=4)
    gas = _array(result["gas_ugm3"], ndim=2)
    bounds = _array(result["bounds_um"], ndim=1)
    if not len(time_h) or number.shape[0] != len(time_h) or mass.shape[0] != len(time_h):
        raise ValueError("run time axis does not match numeric results")
    if number.shape[1] != len(bounds) - 1 or number.shape[2] != len(labels):
        raise ValueError("size distribution axes do not match results")
    if mass.shape[1:3] != number.shape[1:] or mass.shape[3] != len(species):
        raise ValueError("mass axes do not match results")
    if gas.shape != (len(time_h), len(result["gases"])):
        raise ValueError("gas axes do not match results")
    if np.any(bounds <= 0) or np.any(np.diff(bounds) <= 0):
        raise ValueError("invalid size-bin boundaries")

    charts = {}
    total_n = number.sum(axis=(1, 2))
    series = [("total", time_h, total_n)] + [
        (f"C{index+1} {label}", time_h, number[:, :, index].sum(axis=1))
        for index, label in enumerate(labels)]
    charts["charts/overview_number_cm3.png"] = _picture(
        "Particle number over time", "Time (h)", "N (cm^-3)", series)

    dry = [(name, time_h, mass[:, :, :, index].sum(axis=(1, 2)))
           for index, name in enumerate(species) if name != "PH2O"
           and np.any(mass[:, :, :, index] > 0)]
    if dry:
        charts["charts/overview_dry_mass_ugm3.png"] = _picture(
            "Dry aerosol mass over time", "Time (h)", "PM dry (ug m^-3)", dry)
    if "PH2O" in species:
        index = species.index("PH2O")
        charts["charts/overview_particle_water_ugm3.png"] = _picture(
            "Particle water over time", "Time (h)", "PH2O (ug m^-3)",
            [("PH2O", time_h, mass[:, :, :, index].sum(axis=(1, 2)))])
    if gas.shape[1]:
        charts["charts/overview_gases_ugm3.png"] = _picture(
            "Gas concentration over time", "Time (h)", "Gas (ug m^-3)",
            [(name, time_h, gas[:, index]) for index, name in enumerate(result["gases"])
             if np.any(gas[:, index] > 1e-6)]) if np.any(gas > 1e-6) else None
    charts = {name: value for name, value in charts.items() if value is not None}

    met = _array(result["met"], ndim=2)
    if met.shape[0] == len(time_h) and met.shape[1] >= 2:
        charts["charts/overview_temperature_c.png"] = _picture(
            "Air temperature", "Time (h)", "Temperature (C)",
            [("T", time_h, met[:, 0] - 273.15)])
        charts["charts/overview_relative_humidity_pct.png"] = _picture(
            "Relative humidity", "Time (h)", "RH (%)", [("RH", time_h, met[:, 1] * 100)])

    midpoint = np.sqrt(bounds[:-1] * bounds[1:])
    dlog = np.log10(bounds[1:] / bounds[:-1])
    size_series = [("final total", midpoint, number[-1].sum(axis=1) / dlog),
                   ("initial total", midpoint, number[0].sum(axis=1) / dlog)]
    size_series.extend((f"C{index+1} {label}", midpoint, number[-1, :, index] / dlog)
                       for index, label in enumerate(labels))
    final_hour = float(time_h[-1])
    charts["charts/size_number_dNdlogD_cm3.png"] = _picture(
        f"Number size distribution: initial and final t={final_hour:g} h", "Diameter (um)",
        "dN/dlog10(D) (cm^-3)", size_series, logx=True)
    size_mass = [(name, midpoint, mass[-1, :, :, index].sum(axis=1) / dlog)
                 for index, name in enumerate(species) if name != "PH2O"
                 and np.any(mass[-1, :, :, index] > 0)]
    if size_mass:
        charts["charts/size_dry_mass_dMdlogD_ugm3.png"] = _picture(
            f"Dry-mass size distribution: final t={final_hour:g} h", "Diameter (um)",
            "dM/dlog10(D) (ug m^-3)", size_mass, logx=True)

    class_n = number.sum(axis=1)
    total_by_time = class_n.sum(axis=1)
    denominator = np.where(total_by_time != 0, total_by_time, 1.0)
    charts["charts/mix_number_fraction.png"] = _picture(
        "Number fraction by mixing class", "Time (h)", "Fraction",
        [(f"C{index+1} {label}", time_h, class_n[:, index] / denominator)
         for index, label in enumerate(labels)])
    composition = [(name, None, mass[-1, :, :, index].sum(axis=0))
                   for index, name in enumerate(species) if name != "PH2O"
                   and np.any(mass[-1, :, :, index] > 0)]
    if composition:
        charts["charts/mix_dry_mass_by_class_ugm3.png"] = _picture(
            f"Dry mass by class and species: final t={final_hour:g} h", "Mixing class", "Mass (ug m^-3)",
            composition, stacked=True, bar_labels=[f"C{i+1}" for i in range(len(labels))])

    heterogeneity = result.get("het")
    if heterogeneity is not None:
        het = _array(heterogeneity, ndim=2)
        if het.shape[0] == len(time_h) and het.shape[1] >= 7:
            charts["charts/heterogeneous_products_ugm3.png"] = _picture(
                "Heterogeneous products", "Time (h)", "Cumulative (ug m^-3)",
                [("N2O5 uptake", time_h, het[:, 3]), ("HNO3 produced", time_h, het[:, 4]),
                 ("ClNO2 produced", time_h, het[:, 5])])
            charts["charts/heterogeneous_rate_s1.png"] = _picture(
                "N2O5 heterogeneous rate", "Time (h)", "k (s^-1)",
                [("k(N2O5)", time_h, het[:, 0])])
    return charts


def _post_charts(result):
    post = result.get("post") or {}
    charts = {}
    optics = post.get("optics")
    if optics and optics.get("resolved"):
        time_h = _array(optics["time_s"], ndim=1) / 3600.0
        resolved = optics["resolved"]
        internal = optics["internal_mean"]
        charts["charts/optics_extinction_550nm_Mm1.png"] = _picture(
            "550 nm extinction", "Time (h)", "Extinction (Mm^-1)",
            [("resolved", time_h, [row["ext550_Mm"] for row in resolved]),
             ("internal mean", time_h, [row["ext550_Mm"] for row in internal])])
        charts["charts/optics_ssa_550nm.png"] = _picture(
            "550 nm single-scattering albedo", "Time (h)", "SSA",
            [("resolved", time_h, [row["ssa550"] for row in resolved]),
             ("internal mean", time_h, [row["ssa550"] for row in internal])])
        last = resolved[-1]
        charts["charts/optics_spectrum_final_Mm1.png"] = _picture(
            "Final shortwave extinction", "Wavelength (um)", "Extinction (Mm^-1)",
            [("resolved", last["wavelength_sw_um"], last["ext_sw_Mm"]),
             ("internal mean", internal[-1]["wavelength_sw_um"], internal[-1]["ext_sw_Mm"])],
            logx=True)
        if all("bc_abs_enhancement" in row for row in resolved + internal):
            charts["charts/optics_bc_absorption_enhancement.png"] = _picture(
                "BC absorption enhancement", "Time (h)", "Enhancement factor",
                [("resolved", time_h, [row["bc_abs_enhancement"] for row in resolved]),
                 ("internal mean", time_h, [row["bc_abs_enhancement"] for row in internal])])
        if "ext550_by_class_Mm" in last:
            values = last["ext550_by_class_Mm"]
            charts["charts/optics_extinction_by_class_final_Mm1.png"] = _picture(
                "Final 550 nm extinction by class", "Mixing class", "Extinction (Mm^-1)",
                [("extinction", None, values)],
                bar_labels=[f"C{i+1}" for i in range(len(values))])
    ccn = post.get("ccn")
    if ccn and ccn.get("resolved"):
        last = ccn["resolved"][-1]
        mean = ccn["internal_mean"][-1]
        charts["charts/ccn_spectrum_final_cm3.png"] = _picture(
            "Final CCN spectrum", "Supersaturation (%)", "CCN (cm^-3)",
            [("resolved", [row["S_pct"] for row in last["spectrum"]],
              [row["ccn_cm3"] for row in last["spectrum"]]),
             ("internal mean", [row["S_pct"] for row in mean["spectrum"]],
              [row["ccn_cm3"] for row in mean["spectrum"]])], logx=True)
        charts["charts/ccn_droplets_final_cm3.png"] = _picture(
            "Final cloud droplets", "Updraft (m s^-1)", "Nd (cm^-3)",
            [("resolved", [row["w_m_s"] for row in last["activation"]],
              [row["nd_cm3"] for row in last["activation"]]),
             ("internal mean", [row["w_m_s"] for row in mean["activation"]],
              [row["nd_cm3"] for row in mean["activation"]])])
        time_h = _array(ccn["time_s"], ndim=1) / 3600.0
        spectra = last["spectrum"]
        charts["charts/ccn_time_evolution_cm3.png"] = _picture(
            "CCN over time at saved supersaturations", "Time (h)", "CCN (cm^-3)",
            [(f"S={item['S_pct']}%", time_h,
              [row["spectrum"][index]["ccn_cm3"] for row in ccn["resolved"]])
             for index, item in enumerate(spectra)])
        if last.get("activation") and all("actfrac_by_class" in row for row in last["activation"]):
            values = last["activation"]
            charts["charts/ccn_activation_by_class_final.png"] = _picture(
                "Final activated fraction by class", "Mixing class", "Activated fraction",
                [(f"w={item['w_m_s']} m/s", None, item["actfrac_by_class"]) for item in values],
                bar_labels=[f"C{i+1}" for i in range(len(result["class_labels"]))])
        if last.get("kappa_by_class"):
            charts["charts/ccn_kappa_by_class_final.png"] = _picture(
                "Final hygroscopicity by class", "Mixing class", "kappa",
                [("kappa", None, [value or 0 for value in last["kappa_by_class"]])],
                bar_labels=[f"C{i+1}" for i in range(len(result["class_labels"]))])
    hetchem = post.get("hetchem")
    if hetchem and hetchem.get("states"):
        time_h = _array(hetchem["time_s"], ndim=1) / 3600.0
        states = hetchem["states"]
        for gas, rule_modes in states[0]["rules"].items():
            series = []
            for rule in rule_modes["resolved"]:
                for mode, caption in (("resolved", "resolved"), ("internal_mean", "internal mean")):
                    series.append((f"{rule} · {caption}", time_h,
                                   [state["rules"][gas][mode][rule]["k_s1"] for state in states]))
            safe_gas = "".join(ch if ch.isalnum() else "_" for ch in gas)[:32]
            charts[f"charts/hetchem_rate_{safe_gas}_s1.png"] = _picture(
                f"Heterogeneous rate for {gas}", "Time (h)", "k (s^-1)", series)
    return charts


_PARTMC_CHART_NAMES = {
    "charts/partmc_mixing_state.png", "charts/partmc_ccn_spectrum.png",
    "charts/partmc_bc_number.png", "charts/partmc_bc_coating_ratio.png",
}


def _partmc_charts(result):
    truth = result.get("truth") or {}
    if not isinstance(truth, dict) or not truth.get("time_s"):
        return {}
    charts = {}
    hours = _array(truth["time_s"], ndim=1) / 3600.0
    if truth.get("chi"):
        charts["charts/partmc_mixing_state.png"] = _picture(
            "PartMC particle-resolved mixing state", "Time (h)", "Mixing-state index χ",
            [("Particle-resolved χ", hours, truth["chi"])])
    ccn = truth.get("ccn_cm3") or {}
    if ccn:
        supersaturations = sorted((float(key), values) for key, values in ccn.items() if values)
        if supersaturations:
            charts["charts/partmc_ccn_spectrum.png"] = _picture(
                "PartMC particle-resolved CCN at final saved output", "Supersaturation (%)", "CCN (cm⁻³)",
                [("Particle-resolved CCN", [s for s, _ in supersaturations],
                  [values[-1] for _, values in supersaturations])], logx=True)
    if truth.get("n_bc_cm3"):
        charts["charts/partmc_bc_number.png"] = _picture(
            "PartMC BC-containing particle number", "Time (h)", "N_BC (cm⁻³)",
            [("N_BC", hours, truth["n_bc_cm3"])])
    ratios = truth.get("bc_coating_ratio") or []
    valid = [(hour, value) for hour, value in zip(hours, ratios) if value is not None]
    if valid:
        charts["charts/partmc_bc_coating_ratio.png"] = _picture(
            "PartMC BC particle coating ratio", "Time (h)", "Non-BC dry mass / BC mass",
            [("Coating ratio", [point[0] for point in valid], [point[1] for point in valid])])
    return charts


def _timeseries_csv(result):
    out = io.StringIO(newline="")
    writer = csv.writer(out)
    gases = list(result["gases"])
    writer.writerow(["time_s", "time_h", "total_number_cm3", "total_dry_mass_ugm3",
                     "temperature_C", "RH_percent", *[f"gas_{name}_ugm3" for name in gases]])
    number = _array(result["number_m3"], ndim=3) / 1e6
    mass = _array(result["mass_ugm3"], ndim=4)
    species = list(result["species"])
    dry_indices = [i for i, name in enumerate(species) if name != "PH2O"]
    met = _array(result["met"], ndim=2)
    gas = _array(result["gas_ugm3"], ndim=2)
    for index, seconds in enumerate(result["time_s"]):
        writer.writerow([seconds, seconds / 3600, number[index].sum(),
                         mass[index, :, :, dry_indices].sum() if dry_indices else 0,
                         met[index, 0] - 273.15, met[index, 1] * 100,
                         *gas[index].tolist()])
    return out.getvalue().encode("utf-8")


def result_archive(run_id: str, case: dict, result: dict, raw_files: dict[str, bytes]) -> bytes:
    """Return complete JSON plus actual PNG charts in a bounded ZIP archive."""
    raw = _json_bytes(result)
    if len(raw) > MAX_RESULTS_JSON:
        raise ValueError("result is too large for one export; narrow the simulation")
    partmc_charts = _partmc_charts(result) if result.get("model") == "partmc" else {}
    charts = {**_base_charts(result), **_post_charts(result), **partmc_charts}
    manifest = {
        "format": "scrambox-results-1", "run_id": run_id,
        "model": result.get("model", "scram"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "charts": sorted(charts),
        "raw_files": sorted("raw/" + name for name in raw_files),
        "partmc_charts_not_included": sorted(_PARTMC_CHART_NAMES - set(partmc_charts))
        if result.get("model") == "partmc" else [],
        "final_output_time_s": float(_array(result["time_s"], ndim=1)[-1]),
        "final_output_time_h": float(_array(result["time_s"], ndim=1)[-1] / 3600.0),
        "postprocessing_available": sorted((result.get("post") or {}).keys()),
        "postprocessing_not_included": [name for name in ("optics", "ccn", "hetchem")
                                        if not (result.get("post") or {}).get(name)],
        "chart_note": "PNG figures are regenerated from results.json numeric arrays. Size and composition figures use the final saved output, not the web page's current time slider. Missing post-processing is listed and is never computed during export.",
    }
    readme = (
        "SCRAM / PartMC Box result export\n"
        "=======================\n"
        "case.json: submitted case; results.json: UI numeric results rounded to five significant digits; summary.json: run summary.\n"
        "raw/results.npz: original full-precision simulation arrays; raw/grid.json: full grid metadata.\n"
        "raw/case.resolved.json: compiled case; raw/post.json: saved post-processing if present.\n"
        "PartMC runs also include raw/truth.json, raw/partmc_options.json, and raw/out/pmc_0001_*.nc when saved.\n"
        "PartMC particle heatmaps are represented by raw NetCDF; unavailable particle charts are listed in manifest.json.\n"
        "timeseries.csv: selected result time series with explicit units.\n"
        "charts/*.png: regenerated from these numeric results using the units in each axis label.\n"
        f"Final saved output time: {manifest['final_output_time_h']:g} h.\n"
        "Size distributions and class composition use this final time, regardless of the web page slider.\n"
        "Optional post-processing charts appear only when stored post results exist; export does not recompute them.\n"
        "See manifest.json for included charts and omitted post-processing.\n"
    )
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        archive.writestr("manifest.json", _json_bytes(manifest))
        archive.writestr("case.json", _json_bytes(case))
        archive.writestr("results.json", raw)
        archive.writestr("summary.json", _json_bytes(result["summary"]))
        archive.writestr("timeseries.csv", _timeseries_csv(result))
        for name, data in raw_files.items():
            archive.writestr("raw/" + name, data)
        archive.writestr("README.txt", readme.encode("utf-8"))
        for name, image in charts.items():
            archive.writestr(name, image)
    return output.getvalue()
