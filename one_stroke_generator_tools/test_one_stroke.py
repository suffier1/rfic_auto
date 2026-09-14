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

from analyze_sparams import MIXED_MODE, analyze, read_touchstone, target_metrics, contiguous_band
from generate_one_stroke_gds import make_valid_sample, manifest_row, sample_layout, validate_route, route_cells
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
    def test_freeform_size_reproducibility_and_diversity(self):
        from path_checks import simple_path_ok
        for size in (200, 300):
            samples = [make_valid_sample(i, 43, family="serpentine", serpentine_envelope="freeform", size_um=size) for i in range(60)]
            self.assertEqual(len({fingerprint(s) for s in samples}), 60)
            self.assertEqual(Counter(s.mode for s in samples), {"aligned": 20, "offset": 20, "independent": 20})
            self.assertEqual({manifest_row(s)["pri_sweeps"] for s in samples}, {4, 6, 8})
            for s in samples:
                layout = sample_layout(s)
                for route, frame in ((s.pri_route, s.pri_frame), (s.sec_route, s.sec_frame)):
                    self.assertTrue(min(y for _, y in route) <= frame.y_bottom - 20 or max(y for _, y in route) >= frame.y_top + 20)
                for net, layer, side in ((1, 0, "IN"), (2, 1, "OUT")):
                    self.assertTrue(simple_path_ok(layout.cells[net][layer], layout.ports[side + "_P"], layout.ports[side + "_N"]))
                self.assertTrue(all(0 <= x < size / 5 and 0 <= y < size / 5 for layer in (0, 1) for x, y in layout.all_metal(layer)))
            for i in (0, 19, 59):
                self.assertEqual(fingerprint(samples[i]), fingerprint(make_valid_sample(i, 43, family="serpentine", serpentine_envelope="freeform", size_um=size)))
            self.assertNotEqual(fingerprint(samples[0]), fingerprint(make_valid_sample(0, 44, family="serpentine", serpentine_envelope="freeform", size_um=size)))

    def test_freeform_cli_png_and_defect_detection(self):
        import gdstk
        from PIL import Image
        from render_gds import preview_matches, LEFT, TOP
        with tempfile.TemporaryDirectory(prefix="rfic-200-test-") as tmp:
            out = Path(tmp) / "batch"
            run = subprocess.run([sys.executable, str(ROOT / "easy_generate_one_stroke.py"),
                                  "--n", "9", "--seed", "43", "--family", "serpentine", "--size-um", "200",
                                  "--serpentine-envelope", "freeform", "--png", "--outdir", str(out)], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertTrue(json.loads((out / "verification.json").read_text())["png_checked"])
            p = out / "difftx_00000.gds"
            self.assertEqual(verify_one(p, 200, True), [])
            self.assertTrue(verify_one(p, 300, True))
            self.assertTrue(preview_matches(p, p.with_suffix(".png"), 200))
            with Image.open(p.with_suffix(".png")) as image:
                image.putpixel((LEFT, TOP), (0, 0, 0))
                image.save(p.with_suffix(".png"))
            self.assertFalse(preview_matches(p, p.with_suffix(".png"), 200))
            lib = gdstk.read_gds(p)
            cut = next(poly for poly in lib.top_level()[0].polygons if poly.layer == 58)
            cut.translate(.01, 0)
            lib.write_gds(p)
            self.assertTrue(any("VIA cut" in error for error in verify_one(p, 200, True)))

    def test_simple_path_rejects_branch_cycle_and_diagonal(self):
        from path_checks import simple_path_ok
        cells = {(x, 0) for x in range(6)}
        self.assertTrue(simple_path_ok(cells, {(0, 0)}, {(5, 0)}))
        self.assertFalse(simple_path_ok(cells | {(2, 1)}, {(0, 0)}, {(5, 0)}))
        self.assertFalse(simple_path_ok(cells | {(x, 2) for x in range(2, 5)} | {(2, 1), (4, 1)}, {(0, 0)}, {(5, 0)}))
        self.assertFalse(simple_path_ok(cells | {(6, 1)}, {(0, 0)}, {(6, 1)}))

    def test_expanded_envelope_keeps_ports_and_simple_paths(self):
        for index in range(60):
            sample = make_valid_sample(index, 43, family="serpentine", serpentine_envelope="expanded")
            for route, frame, side_x in [(sample.pri_route, sample.pri_frame, 5), (sample.sec_route, sample.sec_frame, 295)]:
                self.assertEqual(route[0], (side_x, frame.y_bottom))
                self.assertEqual(route[-1], (side_x, frame.y_top))
                self.assertTrue(min(y for _, y in route) <= frame.y_bottom - 15 or max(y for _, y in route) >= frame.y_top + 15)
                self.assertGreaterEqual(min(y for _, y in route), 15)
                self.assertLessEqual(max(y for _, y in route), 285)
                validate_route(route)
                cells = route_cells(route)
                ends = {tuple(round(v / 5) for v in p) for p in (route[0], route[-1])}
                for x, y in cells:
                    degree = sum((x + dx, y + dy) in cells for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)])
                    self.assertEqual(degree, 1 if (x, y) in ends else 2)

    def test_mixed_reproducibility_and_nonserpentine_unchanged(self):
        whole = [make_valid_sample(i, 43, family="serpentine", serpentine_envelope="mixed") for i in range(18)]
        self.assertEqual({s.serpentine_envelope for s in whole}, {"bounded", "expanded"})
        split = [make_valid_sample(i, 43, family="serpentine", serpentine_envelope="mixed") for indices in (range(6), range(6, 18)) for i in indices]
        self.assertEqual([fingerprint(s) for s in whole], [fingerprint(s) for s in split])
        for family in ("large_rect", "deep_loop"):
            self.assertEqual(fingerprint(make_valid_sample(4, 8, family=family)), fingerprint(make_valid_sample(4, 8, family=family, serpentine_envelope="expanded")))

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
    def test_common_mode_nonpassivity_is_not_hidden_by_good_dd(self):
        u = np.array([[1, -1, 0, 0], [0, 0, 1, -1], [1, 1, 0, 0], [0, 0, 1, 1]]) / np.sqrt(2)
        mixed = np.zeros((2, 4, 4))
        mixed[:, :2, :2] = [[.1, .8], [.8, .1]]
        mixed[:, 2:, 2:] = np.eye(2) * 1.5
        row = target_metrics(np.array([76e9, 78e9]), u.T @ mixed @ u, 77e9)
        self.assertTrue(row["matching_pass_both_brackets"])
        self.assertFalse(row["full_s_target_passivity_ok"])
        self.assertFalse(row["full_s_sweep_passivity_ok"])
        self.assertFalse(row["matching_pass_both_brackets_with_local_passivity"])
        self.assertAlmostEqual(row["full_s_sigma_max_target"], 1.5)

    def test_band_uses_intersection_not_union(self):
        freq = np.arange(5) * 1e9
        margins = np.array([[-1, -3], [1, -1], [3, 1], [1, 3], [-1, 1]])
        b = contiguous_band(freq, margins, 2, 2)
        self.assertEqual(b["low_est_ghz"], 1.5)
        self.assertEqual(b["high_est_ghz"], 3.5)
        self.assertEqual(b["bw_est_ghz"], 2)
        self.assertFalse(b["left_censored"] or b["right_censored"])

    def test_band_does_not_join_disconnected_passbands(self):
        freq = np.arange(6) * 1e9
        margins = np.array([[1], [1], [-1], [-1], [1], [1]])
        b = contiguous_band(freq, margins, 4, 4)
        self.assertEqual(b["low_est_ghz"], 3.5)
        self.assertTrue(b["right_censored"])
        self.assertFalse(contiguous_band(freq, margins, 1, 2)["supported"])

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
