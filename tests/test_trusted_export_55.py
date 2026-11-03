"""Isolated synthetic completed run; no scientific kernel or production userdirs."""
from __future__ import annotations

import hashlib
import hmac
import io
import json
import os
import tempfile
import threading
import time
import urllib.error
import urllib.request
import zipfile
from http.server import ThreadingHTTPServer
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from scrambox import server


def put_run(workspace, run_id, *, state="done"):
    run = workspace.runs / run_id
    run.mkdir()
    case = {"case": {"name": "synthetic export fixture"}, "kernel": {"isorropia": "fixture"}}
    grid = {"bounds_um": [0.01, 0.1, 1.0], "class_labels": ["A", "B"],
            "species": {"aerosol": ["SO4", "PH2O"], "gases": ["O3"]}, "notes": []}
    values = np.array([[[2123456.789012, 1e6], [1e6, 0.5e6]], [[3e6, 1.5e6], [2e6, 1e6]]])
    mass = np.zeros((2, 2, 2, 2))
    mass[..., 0] = np.array([[[1.0, 0.7], [0.4, 0.2]], [[1.4, 0.9], [0.8, 0.3]]])
    mass[..., 1] = 0.2
    (run / "case.json").write_text(json.dumps(case))
    (run / "case.resolved.json").write_text(json.dumps(case))
    (run / "grid.json").write_text(json.dumps(grid))
    (run / "summary.json").write_text(json.dumps({"cpu_s": 0.0, "n_output": 2, "synthetic": True}))
    (run / "status.json").write_text(json.dumps({"state": state, "name": "synthetic"}))
    np.savez_compressed(run / "results.npz", time_s=np.array([0.0, 3600.0]),
        met=np.array([[298.15, 0.5], [299.15, 0.7]]), number_m3=values,
        mass_ugm3=mass, gas_ugm3=np.array([[1.0], [2.0]]),
        het=np.array([[0.1, 0, 0.2, 1, 2, 3, 0.5], [0.2, 0, 0.3, 2, 4, 6, 0.8]]))
    optics_rows = [{"ext550_Mm": 5.0, "ssa550": 0.9,
                    "bc_abs_enhancement": 1.2, "ext550_by_class_Mm": [3.0, 2.0],
                    "wavelength_sw_um": [0.3, 0.55, 1.0], "ext_sw_Mm": [7.0, 5.0, 2.0]}] * 2
    ccn_rows = [{"spectrum": [{"S_pct": 0.1, "ccn_cm3": 10.0}, {"S_pct": 0.5, "ccn_cm3": 40.0}],
                 "activation": [{"w_m_s": 0.2, "nd_cm3": 11.0, "actfrac_by_class": [0.1, 0.2]},
                                {"w_m_s": 1.0, "nd_cm3": 45.0, "actfrac_by_class": [0.4, 0.5]}],
                 "kappa_by_class": [0.3, 0.4]}] * 2
    het_state = {"rules": {"N2O5": {"resolved": {"dust": {"k_s1": 0.01}},
                                      "internal_mean": {"dust": {"k_s1": 0.02}}}}}
    (run / "post.json").write_text(json.dumps({"optics": {"time_s": [0, 3600],
        "resolved": optics_rows, "internal_mean": optics_rows},
        "ccn": {"time_s": [0, 3600], "resolved": ccn_rows, "internal_mean": ccn_rows},
        "hetchem": {"time_s": [0, 3600], "states": [het_state, het_state]}}))
    return run


def signed_request(port, path, user, secret, *, method="GET", body=b"", cookie=None,
                   skew=0, bad_signature=False):
    timestamp = str(int(time.time()) + skew)
    message = "\n".join((user, timestamp, method, path, hashlib.sha256(body).hexdigest())).encode()
    signature = hmac.new(secret, message, hashlib.sha256).hexdigest()
    headers = {"X-Scram-User": user, "X-Scram-Timestamp": timestamp,
               "X-Scram-Signature": ("0" * 64 if bad_signature else signature)}
    if cookie is not None:
        headers["Cookie"] = cookie
    request = urllib.request.Request(f"http://127.0.0.1:{port}{path}", data=body if method == "POST" else None,
                                     headers=headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            return response.status, response.headers, response.read()
    except urllib.error.HTTPError as error:
        return error.code, error.headers, error.read()


def signed_get(port, path, user, secret):
    return signed_request(port, path, user, secret)


with tempfile.TemporaryDirectory(prefix="scram-export-54-") as temporary:
    original = {key: os.environ.get(key) for key in (
        "SCRAMBOX_USER_DIR", "SCRAMBOX_ACCESS_CODE", "SCRAMBOX_TRUSTED_PROXY_SECRET_FILE")}
    try:
        os.environ["SCRAMBOX_USER_DIR"] = temporary
        os.environ["SCRAMBOX_ACCESS_CODE"] = "synthetic-only-code"
        secret = b"synthetic-export-test-secret-32-bytes-minimum"
        secret_path = Path(temporary) / "trusted_proxy_secret"
        secret_path.write_bytes(secret)
        os.environ["SCRAMBOX_TRUSTED_PROXY_SECRET_FILE"] = str(secret_path)
        owner = server.workspace_for("5002:909090901")
        other = server.workspace_for("5002:909090902")
        put_run(owner, "synthetic_done")
        put_run(owner, "synthetic_queued", state="queued")
        partmc = put_run(owner, "synthetic_partmc_done")
        (partmc / "pmc.spec").write_text("synthetic PartMC fixture\n")
        (partmc / "status.json").write_text(json.dumps({"state": "done", "model": "partmc"}))
        (partmc / "truth.json").write_text(json.dumps({
            "time_s": [0, 3600], "chi": [0.25, 0.45],
            "ccn_cm3": {"0.1": [10, 12], "0.5": [30, 36]},
            "n_bc_cm3": [3, 4], "bc_coating_ratio": [1.1, 1.4],
            "source": "synthetic"}))
        (partmc / "partmc_options.json").write_text(json.dumps({"n_part": 2, "seed": 7}))
        output = partmc / "out"
        output.mkdir()
        with Dataset(output / "pmc_0001_0001.nc", "w", format="NETCDF3_CLASSIC") as netcdf:
            netcdf.createDimension("particle", 2)
            variable = netcdf.createVariable("dry_mass_ug", "f8", ("particle",))
            variable[:] = [1.25, 2.5]
        (output / "notes.txt").write_text("not a result")
        (output / "pmc_0002_0001.nc").write_bytes(b"unrelated trajectory")
        assert not (other.runs / "synthetic_done").exists()
        http = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        worker = threading.Thread(target=http.serve_forever, daemon=True)
        worker.start()
        try:
            session_status, _, session_body = signed_get(http.server_port, "/api/session",
                                                         "5002:909090901", secret)
            assert session_status == 200 and json.loads(session_body) == {
                "mode": "trusted", "user": "5002:909090901"}
            with urllib.request.urlopen(f"http://127.0.0.1:{http.server_port}/api/session", timeout=20) as unsigned_response:
                assert json.loads(unsigned_response.read()) == {"mode": "shared", "user": None}
            cookie_status, _, _ = signed_request(http.server_port, "/api/session",
                "5002:909090901", secret, cookie="sb_session=fake")
            forged_status, _, _ = signed_request(http.server_port, "/api/session",
                "5002:909090901", secret, bad_signature=True)
            expired_status, _, _ = signed_request(http.server_port, "/api/session",
                "5002:909090901", secret, skew=-120)
            login_status, _, _ = signed_request(http.server_port, "/api/login",
                "5002:909090901", secret, method="POST", body=b"{}")
            logout_status, _, _ = signed_request(http.server_port, "/api/logout",
                "5002:909090901", secret, method="POST", body=b"{}")
            assert [cookie_status, forged_status, expired_status, login_status, logout_status] == [403] * 5
            path = "/api/runs/synthetic_done/export"
            status, headers, payload = signed_get(http.server_port, path, "5002:909090901", secret)
            assert status == 200 and headers.get_content_type() == "application/zip"
            assert headers["Content-Disposition"] == 'attachment; filename="scrambox-synthetic_done.zip"'
            with zipfile.ZipFile(io.BytesIO(payload)) as bundle:
                names = set(bundle.namelist())
                assert {"case.json", "results.json", "summary.json", "manifest.json", "README.txt", "timeseries.csv", "raw/results.npz", "raw/grid.json", "raw/case.resolved.json", "raw/post.json"} <= names
                charts = sorted(name for name in names if name.startswith("charts/") and name.endswith(".png"))
                assert len(charts) >= 18, charts
                assert all(bundle.read(name).startswith(b"\x89PNG\r\n\x1a\n") for name in charts)
                results = json.loads(bundle.read("results.json"))
                manifest = json.loads(bundle.read("manifest.json"))
                assert results["time_s"] == [0.0, 3600.0]
                assert results["number_m3"][0][0][0] != 2123456.789012
                raw_npz = np.load(io.BytesIO(bundle.read("raw/results.npz")))
                assert raw_npz["number_m3"][0, 0, 0] == 2123456.789012
                assert manifest["final_output_time_s"] == 3600.0
                assert manifest["final_output_time_h"] == 1.0
                assert not manifest["postprocessing_not_included"]
                assert "final" in bundle.read("README.txt").decode().lower()
                assert b"total_number_cm3" in bundle.read("timeseries.csv")
            foreign_status, _, _ = signed_get(http.server_port, path, "5002:909090902", secret)
            partmc_path = "/api/runs/synthetic_partmc_done/export"
            partmc_status, partmc_headers, partmc_payload = signed_get(
                http.server_port, partmc_path, "5002:909090901", secret)
            assert partmc_status == 200
            assert partmc_headers["Content-Disposition"].endswith('partmcbox-synthetic_partmc_done.zip"')
            with zipfile.ZipFile(io.BytesIO(partmc_payload)) as partmc_zip:
                partmc_names = set(partmc_zip.namelist())
                assert {"raw/truth.json", "raw/partmc_options.json", "raw/out/pmc_0001_0001.nc"} <= partmc_names
                assert "raw/out/pmc_0002_0001.nc" not in partmc_names
                assert "raw/out/notes.txt" not in partmc_names
                assert partmc_zip.read("raw/out/pmc_0001_0001.nc").startswith(b"CDF")
                assert json.loads(partmc_zip.read("manifest.json"))["model"] == "partmc"
                assert not json.loads(partmc_zip.read("manifest.json"))["partmc_charts_not_included"]
                assert len([name for name in partmc_names if name.startswith("charts/partmc_")]) == 4
            foreign_partmc_status, _, _ = signed_get(http.server_port, partmc_path, "5002:909090902", secret)
            foreign_result_status, _, _ = signed_get(
                http.server_port, "/api/runs/synthetic_partmc_done/results", "5002:909090902", secret)
            foreign_case_status, _, _ = signed_get(
                http.server_port, "/api/runs/synthetic_partmc_done/case", "5002:909090902", secret)
            foreign_particles_status, _, _ = signed_get(
                http.server_port, "/api/runs/synthetic_partmc_done/particles", "5002:909090902", secret)
            foreign_delete_status, _, _ = signed_request(
                http.server_port, "/api/runs/synthetic_partmc_done", "5002:909090902", secret, method="DELETE")
            assert {foreign_partmc_status, foreign_result_status, foreign_case_status,
                    foreign_particles_status, foreign_delete_status} == {404}, (foreign_partmc_status, foreign_result_status, foreign_case_status, foreign_particles_status, foreign_delete_status)
            queued_status, _, _ = signed_get(http.server_port, "/api/runs/synthetic_queued/export",
                                             "5002:909090901", secret)
            assert foreign_status == 404 and queued_status == 409
            previous_limit = server.MAX_EXPORT_BODY
            try:
                server.MAX_EXPORT_BODY = 512
                too_large_status, _, _ = signed_get(http.server_port, path, "5002:909090901", secret)
            finally:
                server.MAX_EXPORT_BODY = previous_limit
            assert too_large_status == 413
            try:
                server.run_export(owner, "../bad")
            except KeyError:
                pass
            else:
                raise AssertionError("unsafe run ID was accepted")
            summary = {"owner_http": status, "foreign_http": foreign_status,
                       "trusted_session_http": session_status, "bad_signature_http": forged_status,
                       "signed_cookie_http": cookie_status, "signed_login_http": login_status,
                       "partmc_http": partmc_status, "foreign_partmc_http": foreign_partmc_status,
                       "foreign_particles_http": foreign_particles_status,
                       "queued_http": queued_status, "too_large_http": too_large_status,
                       "charts_png": len(charts),
                       "zip_bytes": len(payload)}
        finally:
            http.shutdown()
            http.server_close()
            worker.join(timeout=5)
    finally:
        for key, value in original.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
assert not Path(temporary).exists()
summary["temporary_root_removed"] = True
print(json.dumps(summary))
