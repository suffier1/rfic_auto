"""두 권선 크기를 일부러 다르게 만든 one-stroke 레이아웃 생성. base 환경 (gdstk).

  python gen_asym.py --n 5000 --seed 2026 --outdir <출력 폴더>

목적: 조건부 모델이 강한 임피던스 변환을 배울 데이터. 기존 데이터는 두 권선 크기가 거의 같아
입출력 임피던스 비가 1:1 근처에 몰려 있다 (DATA_REQUEST.md).

경로는 새로 만들지 않는다. gen_data_colab의 one-stroke 생성기 부품(경로 구성·검증·GDS writer)을
그대로 불러 쓴다. 그 생성기로 만든 1,000장은 협력 랩 EM을 이미 거쳤다
(one-stroke seed 8 corrected 1000장과 verify_asym --baseline으로 대조).
바꾸는 것은 두 권선의 크기를 정하는 곳 하나뿐이다:
  목표 면적비 R = (1차 면적)/(2차 면적)을 [1/3.5, 3.5]에서 로그 균등하게 뽑고,
  두 권선의 폭·높이를 그 비에 맞춘다. 크기는 원래 생성기가 쓰던 범위를 벗어나지 않는다
  (폭 180~250 um, 높이 60~180 um = 포트 피치 규격).

게이트 (하나라도 떨어지면 같은 목표 R로 경로만 다시 뽑는다. 수리는 하지 않는다):
  연결·절연·via·DRC (생성기 기본)
  pin_edge   같은 쪽 두 핀 사이 가장자리 열이 금속으로 이어지지 않는다 (지난 배치 EM 실패 3건의 원인)
  wall       좌우 벽 두 열에 단자가 정확히 2개씩
  단락 턴    한 권선이 제 금속으로 50칸(1,250 um^2) 넘는 면적을 감싸지 않는다 — 권선이 자기
             자신에 닿아 고리가 닫힌 모양. 저장소 공통 정의(checks.self_closed)를 그대로 쓴다.
쓰지 않는 검사 (EM으로 검증된 이 생성기의 1,000장에 돌려 본 결과):
  한붓그리기 경로 검사(path_checks.simple_path_ok) — 352/1000만 통과. 이 생성기의 정상 모양
  (접었다 돌아오는 폭 2칸 손가락, 넓이 없음)을 결함으로 본다.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GEN = HERE.parent / "one_stroke_generator_tools"
sys.path.insert(0, str(GEN)); sys.path.insert(0, str(HERE))

import generate_one_stroke_gds as G                            # noqa: E402
from generation_gates import not_shorted, pin_edge_ok, wall_ok  # noqa: E402,F401

SIZE = 300.0
GRID = 5.0
R_LO, R_HI = 1 / 3.5, 3.5
HEIGHTS = np.arange(60.0, 181.0, 10.0)        # 높이가 10의 배수여야 y_bottom이 5 um 격자에 선다


def width_frame(rng):
    """원래 sample_frame과 같은 폭 범위."""
    xl = G.snap(float(rng.uniform(25.0, 61.0)), GRID)
    xr = G.snap(float(rng.uniform(239.0, 276.0)), GRID)
    return xl, xr


def frames_for_ratio(R, rng, tries=200):
    """목표 면적비 R에 맞는 (1차, 2차) Frame. 2차는 좌우 반전 전 좌표."""
    for _ in range(tries):
        xlp, xrp = width_frame(rng); xls, xrs = width_frame(rng)
        wp, ws = xrp - xlp, xrs - xls
        r = R * ws / wp                                        # 필요한 높이 비 hp/hs
        lo, hi = max(HEIGHTS[0], HEIGHTS[0] / r), min(HEIGHTS[-1], HEIGHTS[-1] / r)
        if lo > hi:
            continue                                            # 이 폭으로는 R을 못 만든다
        hs = float(rng.uniform(lo, hi))
        hs = float(HEIGHTS[np.argmin(abs(HEIGHTS - hs))])
        hp = float(HEIGHTS[np.argmin(abs(HEIGHTS - r * hs))])
        yp, ys = (SIZE - hp) / 2, (SIZE - hs) / 2
        return (G.Frame(xlp, xrp, yp, SIZE - yp), G.Frame(xls, xrs, ys, SIZE - ys))
    return None



def make(index, seed, max_attempts=100):
    """index 하나 -> (Sample, 기록). 목표 R은 index마다 고정하고 경로 난수만 바꿔 재시도한다
    (재시도마다 R을 새로 뽑으면 통과하기 쉬운 R로 분포가 쏠린다)."""
    r0 = np.random.default_rng(np.random.SeedSequence([seed, index, 999_999]))
    R = float(np.exp(r0.uniform(np.log(R_LO), np.log(R_HI))))
    family = G.FAMILIES[index % len(G.FAMILIES)]
    why = {}
    for attempt in range(max_attempts):
        rng = np.random.default_rng(np.random.SeedSequence([seed, index, attempt]))
        fr = frames_for_ratio(R, rng)
        if fr is None:
            why["frame"] = why.get("frame", 0) + 1; continue
        pri_f, sec_f = fr
        try:
            pri = G.build_family_route(family, pri_f, GRID, rng)
            sec = G.mirror_for_right_ports(G.build_family_route(family, sec_f, GRID, rng), SIZE)
            G.validate_route(pri, SIZE); G.validate_route(sec, SIZE)
            s = G.Sample(tag=f"asym_{index:05d}", index=index, seed=seed, attempt=attempt,
                         family=family, mode="asym_ratio", grid_um=GRID, pri_width_um=5.0,
                         sec_width_um=5.0, pri_frame=pri_f, sec_frame=sec_f,
                         pri_route=pri, sec_route=sec, size_um=SIZE)
            lay = G.sample_layout(s)
        except (RuntimeError, ValueError):
            why["route"] = why.get("route", 0) + 1; continue
        for name, ok in (("contract", G.layout_contract_ok(lay)), ("pin_edge", pin_edge_ok(lay)),
                         ("wall", wall_ok(lay)), ("shorted_turn", not_shorted(lay))):
            if not ok:
                why[name] = why.get(name, 0) + 1; break
        else:
            ap = (pri_f.x_right - pri_f.x_left) * (pri_f.y_top - pri_f.y_bottom)
            as_ = (sec_f.x_right - sec_f.x_left) * (sec_f.y_top - sec_f.y_bottom)
            return s, dict(tag=s.tag, index=index, family=family, attempt=attempt,
                           target_ratio=round(R, 4), area_ratio=round(ap / as_, 4),
                           pri_w=pri_f.x_right - pri_f.x_left, pri_h=pri_f.y_top - pri_f.y_bottom,
                           sec_w=sec_f.x_right - sec_f.x_left, sec_h=sec_f.y_top - sec_f.y_bottom,
                           pitch_in_um=pri_f.y_top - pri_f.y_bottom,
                           pitch_out_um=sec_f.y_top - sec_f.y_bottom,
                           route_len_pri=G.route_length(pri), route_len_sec=G.route_length(sec))
    why["target_ratio"] = round(R, 4)
    return None, why


def work(args):
    index, seed, outdir = args
    s, info = make(index, seed)
    if s is None:
        return index, None, info
    G.write_gds(s, Path(outdir) / f"{s.tag}.gds")              # 셀 이름 = 파일 stem = s.tag
    return index, info, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--workers", type=int, default=24)
    a = ap.parse_args()
    gds_dir = Path(a.outdir) / "gds"
    if gds_dir.exists() and any(gds_dir.iterdir()):
        sys.exit(f"{gds_dir}가 비어 있지 않다. 덮어쓰지 않는다.")
    gds_dir.mkdir(parents=True, exist_ok=True)

    rows, fails, idx = [], {}, 0
    with Pool(a.workers) as pool:
        while len(rows) < a.n:
            need = a.n - len(rows)
            batch = [(i, a.seed, str(gds_dir)) for i in range(idx, idx + need)]
            idx += need
            for i, info, why in pool.imap_unordered(work, batch, chunksize=16):
                if info is None:
                    fails[i] = why
                else:
                    rows.append(info)
    rows.sort(key=lambda r: r["index"])
    with open(Path(a.outdir) / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    json.dump(dict(n=len(rows), seed=a.seed, indices_tried=idx, failed_indices=len(fails),
                   ratio_range=[R_LO, R_HI], size_um=SIZE, grid_um=GRID,
                   generator="gen_data_colab/one_stroke_generator_tools/generate_one_stroke_gds.py",
                   fail_reasons={str(k): v for k, v in list(fails.items())[:20]}),
              open(Path(a.outdir) / "gen_meta.json", "w"), indent=1)
    print(f"{gds_dir}  GDS {len(rows)}장  (시도한 index {idx}, 실패한 index {len(fails)})")


if __name__ == "__main__":
    main()
