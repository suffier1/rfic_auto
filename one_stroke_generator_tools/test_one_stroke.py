"""기존 형상 재현, 분할 생성, GDS 규격과 S-parameter 선별의 회귀 검사."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from collections import Counter

import numpy as np

from analyze_sparams import MIXED_MODE, analyze, read_touchstone, target_metrics
from generate_one_stroke_gds import make_valid_sample, manifest_row, sample_layout, validate_route
from verify_one_stroke_gds import verify_one


ROOT = Path(__file__).resolve().parent


def fingerprint(sample):
    layout = sample_layout(sample)
    data = {
        "metal": [sorted(layout.all_metal(k)) for k in [0, 1]],
        "via": sorted(layout.vias),
        "ports": {k: sorted(v) for k, v in sorted(layout.ports.items())},
    }
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


def touchstone(freq, s, form="ri", unit="hz"):
    scale = {"hz": 1, "ghz": 1e9}[unit]
    lines = [f"# {unit} S {form} R 50"]
    for f, matrix in zip(freq, s):
        values = [str(f / scale)]
        for z in matrix.ravel():
            if form == "ri":
                a, b = z.real, z.imag
            else:
                a = abs(z) if form == "ma" else 20 * np.log10(abs(z))
                b = np.angle(z, deg=True)
            values.extend([str(a), str(b)])
        lines.append(" ".join(values))
    return "\n".join(lines)


class GeneratorTests(unittest.TestCase):
    def test_legacy_geometry_fingerprints(self):
        expected = {
            0: "3f42c9224837f7857cb4929dc1cf77c9148ae8d2bc01e5687715840ebb869485",
            1: "5bfe1af58c1e13c05d2ca1839f0cc9293011534f1186f06fd1098fb3f09f219e",
            2: "335dea095c6f9da50243c49c0af52fdd3ef31b78b6d3e4f4bb3922b06cfd6c79",
            8: "ca81c2f7e47e585bc199875545f87e88145c5e1958f5df4b72efb02b532e5ac7",
            26: "23fd20806608ef2c8dc9e274fb9f353e3bd78d13c6b1de4116d2b631ff1d2500",
            999: "dab1ce814e73efe1fdc7a65fafc0651bcbaea5132dc58b8b31d09fcc31daad3f",
        }
        for index, value in expected.items():
            with self.subTest(index=index):
                self.assertEqual(fingerprint(make_valid_sample(index, 8)), value)

    def test_single_family_balances_modes_and_keeps_sweeps(self):
        samples = [make_valid_sample(i, 42, family="serpentine") for i in range(18)]
        self.assertEqual(Counter(s.mode for s in samples), {"aligned": 6, "offset": 6, "independent": 6})
        self.assertEqual({s.family for s in samples}, {"serpentine"})
        for sample in samples:
            row = manifest_row(sample)
            self.assertIn(row["pri_sweeps"], (4, 6))
            self.assertIn(row["sec_sweeps"], (4, 6))
            validate_route(sample.pri_route)
            validate_route(sample.sec_route)

    def test_all_family_schedule(self):
        pairs = Counter((s.family, s.mode) for s in (make_valid_sample(i, 8) for i in range(9)))
        self.assertEqual(len(pairs), 9)
        self.assertEqual(set(pairs.values()), {1})

    def test_cli_chunk_equivalence_and_verification(self):
        with tempfile.TemporaryDirectory(prefix="rfic-test-") as temporary:
            base = Path(temporary)
            manifests = []
            for name, start, count in [("whole", 10000, 6), ("a", 10000, 3), ("b", 10003, 3)]:
                out = base / name
                run = subprocess.run(
                    [sys.executable, str(ROOT / "easy_generate_one_stroke.py"), "--n", str(count),
                     "--seed", "42", "--family", "serpentine", "--start-index", str(start),
                     "--outdir", str(out), "--no-png"], text=True, capture_output=True,
                )
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                self.assertTrue(json.loads((out / "verification.json").read_text())["passed"])
                self.assertEqual(json.loads((out / "dataset_meta.json").read_text())["status"], "complete")
                with (out / "manifest.csv").open() as file:
                    manifests.append(list(csv.DictReader(file)))
                for gds in out.glob("*.gds"):
                    self.assertEqual(verify_one(gds), [])
                self.assertFalse(list(out.glob("*.png")))
            self.assertEqual(manifests[0], manifests[1] + manifests[2])
            for row in manifests[0]:
                name = row["tag"] + ".gds"
                chunk = "a" if int(row["index"]) < 10003 else "b"
                # GDS 내부 timestamp를 제외하고 polygon/label 규격은 검사기로 확인한다.
                import gdstk
                def polygons(path):
                    c = gdstk.read_gds(path).top_level()[0]
                    return sorted((p.layer, p.datatype, tuple(map(tuple, p.points))) for p in c.polygons)
                self.assertEqual(polygons(base / "whole" / name), polygons(base / chunk / name))
            before = (base / "whole" / "manifest.csv").read_bytes()
            collision = subprocess.run(
                [sys.executable, str(ROOT / "easy_generate_one_stroke.py"), "--n", "1", "--outdir", str(base / "whole")],
                capture_output=True,
            )
            self.assertNotEqual(collision.returncode, 0)
            self.assertEqual(before, (base / "whole" / "manifest.csv").read_bytes())

    def test_invalid_cli_args_create_no_output(self):
        with tempfile.TemporaryDirectory(prefix="rfic-invalid-") as temporary:
            out = Path(temporary) / "output"
            for flag, value in [("--n", "0"), ("--seed", "-1"), ("--start-index", "-1"), ("--family", "unknown")]:
                args = [sys.executable, str(ROOT / "generate_one_stroke_gds.py"), "--n", "1", "--seed", "1", "--outdir", str(out), flag, value]
                self.assertNotEqual(subprocess.run(args, capture_output=True).returncode, 0)
                self.assertFalse(out.exists())

    def test_empty_verification_fails_and_saves_report(self):
        with tempfile.TemporaryDirectory(prefix="rfic-empty-") as temporary:
            out = Path(temporary)
            report = out / "verification.json"
            p = subprocess.run([sys.executable, str(ROOT / "verify_one_stroke_gds.py"), str(out), "--report", str(report)], capture_output=True)
            self.assertNotEqual(p.returncode, 0)
            self.assertFalse(json.loads(report.read_text())["passed"])


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.freq = np.array([76e9, 78e9])
        dd = np.array([[[.1, .7+.3j], [.7+.3j, .1]], [[.1, .7-.3j], [.7-.3j, .1]]])
        self.s = np.einsum("pa,fab,bq->fpq", MIXED_MODE.T, dd, MIXED_MODE)

    def test_complex_interpolation_is_not_a_direct_label(self):
        row = target_metrics(self.freq, self.s, 77e9)
        self.assertIsNone(row["matching_pass_direct"])
        self.assertEqual(row["target_source"], "complex_linear_interpolation")
        self.assertTrue(row["matching_pass_both_brackets"])
        self.assertFalse(row["s21_pass_interpolated"])

    def test_exact_target(self):
        row = target_metrics(self.freq, self.s, 76e9)
        self.assertTrue(row["matching_pass_direct"])
        self.assertIsNone(row["matching_pass_interpolated"])
        self.assertEqual(row["left_ghz"], row["right_ghz"])

    def test_no_extrapolation(self):
        for f in [75e9, 79e9]:
            with self.assertRaises(ValueError):
                target_metrics(self.freq, self.s, f)

    def test_touchstone_formats_and_row_order(self):
        # 상호성에 의존하지 않고 비대칭 행렬로 포트의 저장 순서를 확인한다.
        s = np.arange(1, 33).reshape(2, 4, 4) * (.01+.02j)
        for form in ["ri", "ma", "db"]:
            for unit in ["hz", "ghz"]:
                freq, actual = read_touchstone(touchstone(self.freq, s, form, unit))
                np.testing.assert_allclose(freq, self.freq)
                np.testing.assert_allclose(actual, s)

    def test_reject_bad_reference_grid_and_nonfinite(self):
        text = touchstone(self.freq, self.s)
        bad = [text.replace("R 50", "R 75"), text.replace("78000000000.0", "76000000000.0"), text.replace("76000000000.0", "nan"), "[Version] 2.0\n" + text]
        for value in bad:
            with self.assertRaises(ValueError):
                read_touchstone(value)

    def test_missing_and_duplicate_ids(self):
        with tempfile.TemporaryDirectory(prefix="rfic-analysis-") as temporary:
            base = Path(temporary)
            manifest = base / "manifest.csv"
            manifest.write_text("index,tag,family,mode\n0,difftx_00000,serpentine,aligned\n1,difftx_00001,serpentine,offset\n")
            (base / "ostk_00000.s4p").write_text(touchstone(self.freq, self.s))
            rows, report = analyze(manifest, base, 77e9)
            self.assertEqual(report["missing_indices"], [1])
            self.assertEqual(report["all"]["count"], 1)
            self.assertIsNone(report["all"]["matching_pass_direct"])
            (base / "difftx_00000.s4p").write_text(touchstone(self.freq, self.s))
            with self.assertRaises(ValueError):
                analyze(manifest, base, 77e9)


if __name__ == "__main__":
    unittest.main()
