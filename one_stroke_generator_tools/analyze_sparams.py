#!/usr/bin/env python3
"""manifest와 s4p를 숫자 ID로 연결하고 지정 주파수의 차동 성능을 집계한다."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import tarfile
from pathlib import Path

import numpy as np


MIXED_MODE = np.array([[1, -1, 0, 0], [0, 0, 1, -1]]) / np.sqrt(2.0)


def read_touchstone(text: str) -> tuple[np.ndarray, np.ndarray]:
    """Touchstone 1의 4-port S를 읽는다. 포트별 기준 임피던스는 50 Ω이다."""

    options = None
    values = []
    for line in text.splitlines():
        line = line.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            raise ValueError("Touchstone 2는 지원하지 않습니다. Touchstone 1 s4p로 export하세요.")
        if line.startswith("#"):
            if options is not None:
                raise ValueError("중복 Touchstone option 행")
            options = line[1:].lower().split()
        else:
            values.extend(float(v.replace("D", "E").replace("d", "e")) for v in line.split())
    if options is None or len(options) != 5 or options[1] != "s" or options[3] != "r":
        raise ValueError("# hz S ri R 50 형태의 option 행이 필요합니다.")
    if float(options[4]) != 50:
        raise ValueError("이 선별 도구는 포트별 기준 임피던스 50 Ω만 지원합니다.")
    scales = {"hz": 1, "khz": 1e3, "mhz": 1e6, "ghz": 1e9}
    if options[0] not in scales or options[2] not in ("ri", "ma", "db"):
        raise ValueError(f"지원하지 않는 Touchstone 옵션: {options}")
    if not values or len(values) % 33:
        raise ValueError("4-port 데이터는 주파수당 33개의 숫자가 필요합니다.")
    data = np.asarray(values).reshape(-1, 33)
    freq = data[:, 0] * scales[options[0]]
    a, b = data[:, 1::2], data[:, 2::2]
    if options[2] == "ri":
        s = a + 1j * b
    else:
        magnitude = a if options[2] == "ma" else 10.0 ** (a / 20.0)
        s = magnitude * np.exp(1j * np.deg2rad(b))
    if not np.isfinite(data).all() or not np.isfinite(s).all() or np.any(np.diff(freq) <= 0):
        raise ValueError("비유한 값 또는 중복/역순 주파수")
    # 4-port Touchstone 1은 S11,S12,...,S21,S22,... 순서이다.
    return freq, s.reshape(-1, 4, 4)


def iter_s4p(path: Path):
    """디렉터리 또는 tar.gz에서 읽기만 수행한다. 압축을 해제하지 않는다."""

    if path.is_dir():
        for file in sorted(path.rglob("*.s4p")):
            yield str(file), file.read_text(encoding="utf-8-sig")
    else:
        with tarfile.open(path, "r|*") as archive:
            for member in archive:
                if member.isfile() and member.name.lower().endswith(".s4p"):
                    yield member.name, archive.extractfile(member).read().decode("utf-8-sig")


def differential_db(s: np.ndarray) -> np.ndarray:
    dd = np.einsum("ap,...pq,bq->...ab", MIXED_MODE, s, MIXED_MODE)
    return 20 * np.log10(np.maximum(np.abs(dd), 1e-300))


def pass_flags(db: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    s21 = db[..., 1, 0] >= -3
    return s21, s21 & (db[..., 0, 0] <= -10) & (db[..., 1, 1] <= -10)


def contiguous_band(freq: np.ndarray, margins: np.ndarray, left: int, right: int) -> dict:
    """타깃을 지지하는 연속 샘플 구간과 문턱 교차점의 선형 dB 추정치를 구한다.

    margins의 모든 열이 0 이상이어야 통과한다. 타깃이 샘플 사이에 있으면
    양옆 지점 모두의 통과를 요구한다. sweep 끝에 닿으면 실제 대역 끝은 미확정이다.
    """

    good = (margins >= 0).all(axis=1)
    result = {"supported": bool(good[left] and good[right]),
              "low_sample_ghz": None, "high_sample_ghz": None, "sample_span_ghz": None,
              "low_est_ghz": None, "high_est_ghz": None, "bw_est_ghz": None,
              "left_censored": None, "right_censored": None}
    if not result["supported"]:
        return result
    lo, hi = left, right
    while lo > 0 and good[lo - 1]:
        lo -= 1
    while hi + 1 < len(freq) and good[hi + 1]:
        hi += 1
    low, high = float(freq[lo]), float(freq[hi])
    if lo > 0:
        crossings = []
        for bad, passed in zip(margins[lo - 1], margins[lo]):
            if bad < 0:
                crossings.append(freq[lo - 1] + (freq[lo] - freq[lo - 1]) * (-bad) / (passed - bad))
        low = float(max(crossings))
    if hi + 1 < len(freq):
        crossings = []
        for passed, bad in zip(margins[hi], margins[hi + 1]):
            if bad < 0:
                crossings.append(freq[hi] + (freq[hi + 1] - freq[hi]) * passed / (passed - bad))
        high = float(min(crossings))
    result.update(low_sample_ghz=float(freq[lo] / 1e9), high_sample_ghz=float(freq[hi] / 1e9),
                  sample_span_ghz=float((freq[hi] - freq[lo]) / 1e9),
                  low_est_ghz=low / 1e9, high_est_ghz=high / 1e9, bw_est_ghz=(high - low) / 1e9,
                  left_censored=lo == 0, right_censored=hi == len(freq) - 1)
    return result


def target_metrics(freq: np.ndarray, s: np.ndarray, target: float) -> dict:
    """직접 EM과 보간을 구분하며, 양옆 지점 통과를 대역 전체 통과로 부르지 않는다."""

    exact = np.flatnonzero(np.isclose(freq, target, rtol=0, atol=1))
    if len(exact):
        left = right = int(exact[0])
        value = s[left]
        source = "direct_em"
    else:
        right = int(np.searchsorted(freq, target))
        if right == 0 or right == len(freq):
            raise ValueError("목표 주파수가 EM 범위 밖에 있습니다. 외삽하지 않습니다.")
        left = right - 1
        weight = (target - freq[left]) / (freq[right] - freq[left])
        value = (1 - weight) * s[left] + weight * s[right]
        source = "complex_linear_interpolation"
    db = differential_db(value)
    all_db = differential_db(s)
    p, q = pass_flags(db)
    bp, bq = pass_flags(all_db[[left, right]])
    peak = int(np.argmax(all_db[:, 1, 0]))
    rl_margins = np.column_stack((-10 - all_db[:, 0, 0], -10 - all_db[:, 1, 1]))
    joint_margins = np.column_stack((rl_margins, all_db[:, 1, 0] + 3))
    # 실수 50 Ω 기준의 수동 4-port는 모든 특이값이 1 이하여야 한다.
    # 차동 S21 통과와 전체 S행렬의 물리적 정상성은 별개로 기록한다.
    sigma = np.linalg.svd(s, compute_uv=False)[:, 0]
    target_sigma = float(np.linalg.svd(value, compute_uv=False)[0])
    brackets_passive = bool((sigma[[left, right]] <= 1 + 1e-6).all())
    row = {
        "target_ghz": target / 1e9,
        "target_source": source,
        "left_ghz": float(freq[left] / 1e9),
        "right_ghz": float(freq[right] / 1e9),
        "sdd21_db": float(db[1, 0]),
        "sdd11_db": float(db[0, 0]),
        "sdd22_db": float(db[1, 1]),
        "s21_pass_direct": bool(p) if len(exact) else None,
        "matching_pass_direct": bool(q) if len(exact) else None,
        "s21_pass_interpolated": bool(p) if not len(exact) else None,
        "matching_pass_interpolated": bool(q) if not len(exact) else None,
        "s21_pass_both_brackets": bool(bp.all()),
        "matching_pass_both_brackets": bool(bq.all()),
        "peak_frequency_ghz": float(freq[peak] / 1e9),
        "peak_at_sweep_edge": peak in (0, len(freq) - 1),
        "target_below_peak_db": float(all_db[peak, 1, 0] - db[1, 0]),
        "full_s_sigma_max_target": target_sigma,
        "full_s_target_passivity_ok": target_sigma <= 1 + 1e-6,
        "full_s_brackets_passivity_ok": brackets_passive,
        "full_s_sweep_passivity_ok": bool((sigma <= 1 + 1e-6).all()),
        "full_s_sweep_sigma_max": float(sigma.max()),
        "full_s_reciprocity_error_max": float(np.abs(s - s.swapaxes(-1, -2)).max()),
        "matching_pass_both_brackets_with_local_passivity": bool(bq.all() and brackets_passive),
    }
    for prefix, margins in (("rl_band", rl_margins), ("joint_band", joint_margins)):
        row.update({f"{prefix}_{key}": val for key, val in contiguous_band(freq, margins, left, right).items()})
    return row


def summarize(rows: list[dict]) -> dict:
    direct = [r for r in rows if r["target_source"] == "direct_em"]
    estimated = [r for r in rows if r["target_source"] != "direct_em"]
    return {
        "count": len(rows),
        "direct_count": len(direct),
        "interpolated_count": len(estimated),
        "s21_pass_direct": sum(r["s21_pass_direct"] for r in direct) if direct else None,
        "matching_pass_direct": sum(r["matching_pass_direct"] for r in direct) if direct else None,
        "s21_pass_interpolated": sum(r["s21_pass_interpolated"] for r in estimated) if estimated else None,
        "matching_pass_interpolated": sum(r["matching_pass_interpolated"] for r in estimated) if estimated else None,
        "s21_pass_both_brackets": sum(r["s21_pass_both_brackets"] for r in rows),
        "matching_pass_both_brackets": sum(r["matching_pass_both_brackets"] for r in rows),
        "full_s_target_nonpassive_count": sum(not r["full_s_target_passivity_ok"] for r in rows),
        "full_s_brackets_nonpassive_count": sum(not r["full_s_brackets_passivity_ok"] for r in rows),
        "full_s_sweep_nonpassive_count": sum(not r["full_s_sweep_passivity_ok"] for r in rows),
        "matching_pass_both_brackets_with_local_passivity": sum(r["matching_pass_both_brackets_with_local_passivity"] for r in rows),
    }


def analyze(manifest_path: Path, sparams: Path, target: float) -> tuple[list[dict], dict]:
    with manifest_path.open(newline="", encoding="utf-8") as file:
        manifest_rows = list(csv.DictReader(file))
    manifest = {int(r["index"]): r for r in manifest_rows}
    if not manifest or len(manifest) != len(manifest_rows):
        raise ValueError("manifest가 비어 있거나 index가 중복됩니다.")
    rows, seen = [], set()
    for name, text in iter_s4p(sparams):
        match = re.search(r"_(\d+)$", Path(name).stem)
        if not match:
            raise ValueError(f"숫자 ID를 읽을 수 없는 이름: {name}")
        index = int(match[1])
        if index not in manifest or index in seen:
            raise ValueError(f"manifest에 없는 ID 또는 중복 s4p: {name}")
        seen.add(index)
        freq, s = read_touchstone(text)
        rows.append({**manifest[index], "s4p": name, **target_metrics(freq, s, target)})
    if not rows:
        raise ValueError("s4p 파일이 없습니다.")
    rows.sort(key=lambda r: int(r["index"]))
    report = {
        "manifest": str(manifest_path.resolve()),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "sparams": str(sparams.resolve()),
        "target_ghz": target / 1e9,
        "reference": "50 ohm per single-ended port; 100 ohm differential; pairs (1,2),(3,4)",
        "criteria": "Sdd21 >= -3 dB; matching also requires Sdd11 and Sdd22 <= -10 dB at the same frequency",
        "note": "Brackets are sampled points. Target S is complex-linearly interpolated when absent. Band edges use linear dB margins; censored edges give only a sweep-limited width. RL band uses both reflections; joint band also requires Sdd21 >= -3 dB.",
        "quality_note": "Raw threshold counts and bandwidths do not certify EM validity. Full real-50-ohm S passivity requires sigma_max <= 1 + 1e-6; inspect target, brackets and full sweep separately. No passivity enforcement or data correction is applied. Passing this necessary test does not validate calibration or mesh convergence.",
        "expected_count": len(manifest),
        "missing_indices": sorted(set(manifest) - seen),
        "bracket_pairs_ghz": sorted({(r["left_ghz"], r["right_ghz"]) for r in rows}),
        "all": summarize(rows),
        "by_family": {f: summarize([r for r in rows if r["family"] == f]) for f in sorted({r["family"] for r in rows})},
        "by_family_mode": {
            f"{f}/{m}": summarize([r for r in rows if r["family"] == f and r["mode"] == m])
            for f, m in sorted({(r["family"], r["mode"]) for r in rows})
        },
    }
    return rows, report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--sparams", type=Path, required=True, help="s4p 디렉터리 또는 tar.gz")
    parser.add_argument("--freq-ghz", type=float, default=77)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if not np.isfinite(args.freq_ghz) or args.freq_ghz <= 0:
        parser.error("freq-ghz는 양의 유한한 값이어야 합니다.")
    if args.outdir.exists():
        raise FileExistsError(f"기존 분석 결과를 덮어쓰지 않습니다: {args.outdir}")
    rows, report = analyze(args.manifest, args.sparams, args.freq_ghz * 1e9)
    args.outdir.mkdir(parents=True)
    with (args.outdir / "per_sample.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.outdir / "summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
