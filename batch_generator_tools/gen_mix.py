"""mix 배치: 크기가 비슷한 P/S + 경로 다양화 3,000장, 비대칭 2,000장. base 환경 (gdstk).

  python gen_mix.py --seed 2027 --outdir <출력 폴더>
  python gen_mix.py --pilot --seed 2027 --outdir <출력 폴더>   # 구성별 대표 몇 장

구성 (index로 결정론적 배정, 재시도해도 배정은 유지):
  index 0..2999   similar  family = index%3, 크기급 S/M/L = (index//3)%3, 상대배치 aligned/offset/free = (index//9)%3
                  P/S 폭·높이 각각 min/max >= 0.8 — frame뿐 아니라 실제 금속(각 권선이 제 벽에서 뻗은 길이,
                  y 범위)에도 건다. 크기급은 P frame 면적의 로그 3등분.
                    aligned: 레이아웃에서 S 외곽 = P 외곽 (dx = 0). 경로 dial은 따로 뽑는다 (serpentine은 행 수만
                             같고 행 간격은 따로 뽑아 겹침을 강제하지 않는다)
                    offset : S 치수 = P 치수, 레이아웃 중심을 x로 10 um 이상 옮긴다 (가능한 범위에서 균등)
                    free   : S 치수·위치를 0.8 조건 안에서 따로 뽑는다
                  (S frame은 반전 전 좌표로 만들어지므로 원하는 레이아웃 자리의 거울상을 준다)
  index 3000..3999 asym_lo  목표 면적비 R = P/S 로그균등 [1/3.5, 0.8], family = index%3
  index 4000..4999 asym_hi  R 로그균등 [1.25, 3.5]
                  폭 비를 먼저 R과 같은 쪽으로 뽑고 나머지를 높이 비로 채운다 (폭 범위 180~250 um가 좁아
                  큰 비는 결국 높이에서 나온다). snap 뒤 frame 면적비가 구간 밖이면 다시 뽑는다.

경로: one-stroke 생성기의 부품(경로 검증, 레이아웃 조립, 계약 게이트, GDS writer)을 그대로 쓰되,
경로 구성은 dial을 밖으로 뺀 새 builder를 쓴다. 굴곡 위치는 서로 2칸 이상 떨어뜨려 폭 2칸으로 접히는
구간(finger)과 3x3 금속 블록이 생기지 않게 한다.
  열린 루프(large_rect/deep_loop): 아래·위·오른쪽 변의 굴곡 깊이(격자 단위)와 굴곡 수를 직접 뽑는다.
      large_rect 깊이 0~6칸(높이가 허용하는 만큼), deep_loop 3칸~(높이/2-3칸). 띠 사이에 빈 행을 항상 남긴다.
  serpentine: 왕복 4/6/8행(평균 행 간격 3칸 이상이 되는 수만), 행 간격 불균일(각 >=2칸), 행마다 반환 위치
      (폭의 30%까지), 굴곡은 min(2, 행 간격-2)칸 이하로 제한해 다음 행에 닿지 않는다. 마지막 행은 안쪽으로 굽힌다.

게이트 (하나라도 떨어지면 버리고 같은 배정으로 다시 뽑는다. 수리 없음):
  연결·절연·via·DRC (생성기 계약), pin_edge, wall(벽 단자 2개씩),
  no_hole: 어느 권선도 제 금속으로 구멍을 감싸지 않는다 — 크기 무관 (기존 50칸 기준보다 엄격).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
GEN = ROOT / "one_stroke_generator_tools"
sys.path.insert(0, str(GEN)); sys.path.insert(0, str(HERE))

import generate_one_stroke_gds as G                                        # noqa: E402
from em_gds_contract import M8, M9, N, NET_IN, NET_OUT                      # noqa: E402,F401
from gen_asym import GRID, HEIGHTS, SIZE                                    # noqa: E402
from generation_gates import hole_cells, no_hole, pin_edge_ok, wall_ok     # noqa: E402,F401

PORT_X = 5.0
W_MIN, W_MAX = 180.0, 250.0                    # 원래 생성기의 x 범위에서 나오는 폭 (snap 후)
X_LEFT = (25.0, 61.0); X_RIGHT = (239.0, 276.0)
SIZE_EDGES = np.geomspace(W_MIN * HEIGHTS[0], W_MAX * HEIGHTS[-1], 4)      # P 면적 로그 3등분
SIZE_CLASSES = ("S", "M", "L")
RELS = ("aligned", "offset", "free")
GROUPS = {"similar": None, "asym_lo": (1 / 3.5, 0.8), "asym_hi": (1.25, 3.5)}
N_SIMILAR, N_ASYM = 3000, 1000


# ----------------------------------------------------------------- 배정
def assignment(index: int) -> dict:
    fam = G.FAMILIES[index % 3]
    if index < N_SIMILAR:
        return dict(group="similar", family=fam, size_class=SIZE_CLASSES[(index // 3) % 3],
                    rel=RELS[(index // 9) % 3])
    grp = "asym_lo" if index < N_SIMILAR + N_ASYM else "asym_hi"
    return dict(group=grp, family=fam, size_class="", rel="free")


def pilot_indices(per_similar=2, per_asym=4):
    """구성 칸마다 대표 index 몇 개 (similar 27칸, asym 6칸)."""
    out, seen = [], {}
    for i in range(N_SIMILAR + 2 * N_ASYM):
        a = assignment(i); k = (a["group"], a["family"], a["size_class"], a["rel"])
        cap = per_similar if a["group"] == "similar" else per_asym
        if seen.get(k, 0) < cap:
            seen[k] = seen.get(k, 0) + 1; out.append(i)
    return out


# ----------------------------------------------------------------- frame
def dims_ok(a: float, b: float) -> bool:
    return min(a, b) / max(a, b) >= 0.8


def size_class_of(w: float, h: float) -> str:
    i = int(np.clip(np.searchsorted(SIZE_EDGES, w * h, side="right") - 1, 0, 2))
    return SIZE_CLASSES[i]


def centered(xl: float, xr: float, h: float) -> G.Frame:
    return G.Frame(xl, xr, (SIZE - h) / 2, (SIZE + h) / 2)


def p_frame_in_class(cls: str, rng) -> G.Frame:
    lo, hi = SIZE_EDGES[SIZE_CLASSES.index(cls)], SIZE_EDGES[SIZE_CLASSES.index(cls) + 1]
    for _ in range(1000):
        xl = G.snap(float(rng.uniform(*X_LEFT)), GRID); xr = G.snap(float(rng.uniform(*X_RIGHT)), GRID)
        h = float(rng.choice(HEIGHTS)); a = (xr - xl) * h
        if lo <= a < hi or (cls == "L" and a >= hi - 1e-6):
            return centered(xl, xr, h)
    raise RuntimeError("size class frame")


def x_range_for_width(w: float):
    """폭 w인 frame의 x_left가 갈 수 있는 범위 (원래 생성기의 x 범위 안, 5 um 격자 위)."""
    lo, hi = max(X_LEFT[0], X_RIGHT[0] - w), min(X_LEFT[1], X_RIGHT[1] - w)
    return float(np.ceil(lo / GRID) * GRID), float(np.floor(hi / GRID) * GRID)


def similar_partner(pf: G.Frame, rel: str, rng):
    """P와 폭·높이가 각각 0.8 이상 비슷한 S frame (좌우 반전 전 좌표)."""
    wp, hp = pf.x_right - pf.x_left, pf.y_top - pf.y_bottom
    # S는 나중에 좌우 반전되므로, 레이아웃에서 원하는 자리의 거울상을 반전 전 frame으로 준다.
    if rel == "aligned":                     # 레이아웃에서 S 외곽 = P 외곽
        return G.Frame(SIZE - pf.x_right, SIZE - pf.x_left, pf.y_bottom, pf.y_top)
    if rel == "offset":                      # 같은 치수, 레이아웃 중심을 x로 |dx| >= 2칸 옮긴다
        lo, hi = x_range_for_width(wp)       # 반전 전 x_left 범위. 레이아웃 x_left = SIZE - x_left - wp
        cands = [x for x in np.arange(lo, hi + 0.1, GRID)
                 if abs((SIZE - x - wp) - pf.x_left) >= 2 * GRID - 1e-6]
        if not cands:
            return None
        xl = float(rng.choice(cands))
        return centered(xl, xl + wp, hp)
    for _ in range(100):
        ws = G.snap(float(rng.uniform(max(W_MIN, 0.8 * wp), min(W_MAX, 1.25 * wp))), GRID)
        hs = float(rng.choice([h for h in HEIGHTS if dims_ok(h, hp)]))
        lo, hi = x_range_for_width(ws)
        xl = G.snap(float(rng.uniform(lo, hi)), GRID)
        if dims_ok(wp, ws) and lo - 1e-6 <= xl <= hi + 1e-6:
            return centered(xl, xl + ws, hs)
    return None


def asym_frames(R: float, lo_bin: float, hi_bin: float, rng, tries=300):
    """목표 면적비 R(P/S). 폭 비 r_w를 R과 같은 쪽에서 로그균등하게 먼저 뽑고 높이 비로 나머지를 채운다.
    snap 뒤 frame 면적비가 [lo_bin, hi_bin] 밖이면 다시 뽑는다."""
    w_span = np.log(W_MAX / W_MIN)
    for _ in range(tries):
        wp = G.snap(float(rng.uniform(W_MIN, W_MAX)), GRID)
        lim = float(np.clip(np.log(R), -w_span, w_span))          # 폭만으로 낼 수 있는 최대 |log 비|
        lrw = float(rng.uniform(0, lim)) if lim >= 0 else float(rng.uniform(lim, 0))
        ws = G.snap(wp / np.exp(lrw), GRID)
        if not (W_MIN <= ws <= W_MAX):
            continue
        r_h = R * ws / wp                                          # 필요한 높이 비 hp/hs
        lo, hi = max(HEIGHTS[0], HEIGHTS[0] / r_h), min(HEIGHTS[-1], HEIGHTS[-1] / r_h)
        if lo > hi:
            continue
        hs = float(HEIGHTS[np.argmin(abs(HEIGHTS - rng.uniform(lo, hi)))])
        hp = float(HEIGHTS[np.argmin(abs(HEIGHTS - r_h * hs))])
        if not (lo_bin - 1e-9 <= wp * hp / (ws * hs) <= hi_bin + 1e-9):
            continue
        xlp = G.snap(float(rng.uniform(*x_range_for_width(wp))), GRID)
        xls = G.snap(float(rng.uniform(*x_range_for_width(ws))), GRID)
        return centered(xlp, xlp + wp, hp), centered(xls, xls + ws, hs)
    return None


# ----------------------------------------------------------------- 경로 builder
def breaks(x0: float, x1: float, n: int, rng, sep: int = 2):
    """[x0, x1] 안쪽(양 끝 2칸 제외)에서 서로 sep칸 이상 떨어진 굴곡 위치 n개 (부족하면 되는 만큼)."""
    lo, hi = sorted((x0, x1))
    cands = list(np.arange(lo + 2 * GRID, hi - 2 * GRID + 0.1 * GRID, GRID))
    rng.shuffle(cands)
    out = []
    for c in cands:
        if len(out) == n:
            break
        if all(abs(c - o) >= sep * GRID - 1e-6 for o in out):
            out.append(float(c))
    out.sort()
    return out if x1 > x0 else out[::-1]


def band_h(start, end, y_min, y_max, rng, n_breaks):
    """x는 한 방향, y는 [y_min, y_max] 띠 안에서 굴곡마다 자유롭게."""
    x0, y0 = start; x1, _ = end
    levels = np.arange(G.snap(y_min, GRID), G.snap(y_max, GRID) + 0.1 * GRID, GRID)
    pts, cur = [start], y0
    for x in breaks(x0, x1, n_breaks, rng):
        G.append_point(pts, (x, cur)); cur = float(rng.choice(levels)); G.append_point(pts, (x, cur))
    G.append_point(pts, (x1, cur)); G.append_point(pts, end)
    return G.compact_route(pts)


def band_v(start, end, x_min, x_max, rng, n_breaks):
    x0, y0 = start; _, y1 = end
    levels = np.arange(G.snap(x_min, GRID), G.snap(x_max, GRID) + 0.1 * GRID, GRID)
    pts, cur = [start], x0
    for y in breaks(y0, y1, n_breaks, rng):
        G.append_point(pts, (cur, y)); cur = float(rng.choice(levels)); G.append_point(pts, (cur, y))
    G.append_point(pts, (cur, y1)); G.append_point(pts, end)
    return G.compact_route(pts)


def jag_h(start, end, sign: int, rng, max_off: int, n_breaks: int):
    """작은 직교 굴곡 (generator.jagged_horizontal과 같되 굴곡 위치가 2칸 이상 떨어진다)."""
    x0, y0 = start; x1, _ = end
    pts, cur = [start], y0
    for x in breaks(x0, x1, n_breaks, rng):
        G.append_point(pts, (x, cur)); cur = y0 + sign * int(rng.integers(0, max_off + 1)) * GRID
        G.append_point(pts, (x, cur))
    G.append_point(pts, (x1, cur)); G.append_point(pts, end)
    return G.compact_route(pts)


def build_loop(frame: G.Frame, rng, family: str):
    """열린 루프. 아래(좌→우)·오른쪽(아래→위)·위(우→좌) 세 띠가 서로 빈 행/열을 사이에 두므로
    자기 접촉이 구조적으로 없다. 깊이 d_*는 격자 칸 수 (뽑은 한계), *_real은 실제로 쓰인 최대 깊이."""
    xl, xr, yb, yt = frame.x_left, frame.x_right, frame.y_bottom, frame.y_top
    hg, wg = int(round((yt - yb) / GRID)), int(round((xr - xl) / GRID))
    if family == "large_rect":
        hi = max(3, min(6, hg // 2 - 3))                # 0~30 um. deep_loop(3칸 이상)와 3~6칸이 겹친다
        d_b, d_t = int(rng.integers(0, hi + 1)), int(rng.integers(0, hi + 1))
        d_r = int(rng.integers(0, min(6, int(0.45 * wg)) + 1))
        nb = [int(rng.integers(2, 7)) for _ in range(3)]
    else:
        hi = max(3, hg // 2 - 3)
        d_b, d_t = int(rng.integers(3, hi + 1)), int(rng.integers(3, hi + 1))
        d_r = int(rng.integers(3, max(3, int(0.45 * wg)) + 1))
        nb = [int(rng.integers(3, 9)) for _ in range(3)]
    # 오른쪽 띠가 굴곡을 넣을 자리(7행)를 남긴다. 낮은 frame이면 오른쪽은 직선으로 둔다.
    while d_r > 0 and hg - d_b - d_t < 7:
        if max(d_b, d_t) > 3:
            if d_b >= d_t: d_b -= 1
            else: d_t -= 1
        else:
            d_r = 0
    while hg - d_b - d_t < 3:
        d_t -= 1
    y_low, y_high = yb + (d_b + 1) * GRID, yt - (d_t + 2) * GRID
    bottom = band_h((xl, yb), (xr, yb), yb, yb + d_b * GRID, rng, nb[0])
    if d_r > 0:
        right = band_v((xr, y_low), (xr, y_high), xr - d_r * GRID, xr, rng, nb[2])
    else:
        right = [(xr, y_low), (xr, y_high)]
    top = band_h((xr, yt), (xl, yt), yt - d_t * GRID, yt, rng, nb[1])
    route = G.join_routes([(PORT_X, yb), (xl, yb)], bottom, [(xr, yb), (xr, y_low)], right,
                          [(xr, y_high), (xr, yt)], top, [(xl, yt), (PORT_X, yt)])
    real = dict(db_real=int(round((max(y for _, y in bottom) - yb) / GRID)),
                dt_real=int(round((yt - min(y for _, y in top)) / GRID)),
                dr_real=int(round((xr - min(x for x, _ in right)) / GRID)))
    return route, dict(d_b=d_b, d_t=d_t, d_r=d_r, nb=nb, **real)


def build_serp(frame: G.Frame, rng, fixed: dict | None = None):
    """왕복 n행. 행 간격 gaps(칸, 각 >=2, 합 = 높이), 행 i의 반환 위치 ret[i](끝에서 안쪽으로 몇 칸),
    행 i의 굴곡 높이 bump[i] <= gaps[i]-2 (다음 행과 빈 행 하나를 남긴다). 마지막 행은 아래로 굽힌다.
    fixed={'n': k}면 행 수만 맞추고 간격은 따로 뽑는다."""
    xl, xr, yb, yt = frame.x_left, frame.x_right, frame.y_bottom, frame.y_top
    hg, wg = int(round((yt - yb) / GRID)), int(round((xr - xl) / GRID))
    # 행 간격이 전부 최소(2칸)로 고정되면 같은 크기의 샘플끼리 거의 같은 모양이 된다.
    # 평균 간격이 3칸 이상 되는 행 수만 허용한다 (n=6은 높이 80 um부터, n=8은 110 um부터).
    allowed = [k for k in (4, 6, 8) if (k - 1) * 3 <= hg]
    n = fixed["n"] if fixed and fixed["n"] in allowed else int(rng.choice(allowed))
    gaps = (2 + rng.multinomial(hg - 2 * (n - 1), rng.dirichlet(np.ones(n - 1)))).tolist()
    K = max(0, int(0.3 * wg))
    ret = [int(rng.integers(0, K + 1)) for _ in range(n - 1)]
    bump = [int(rng.integers(0, min(2, g - 2) + 1)) for g in gaps]
    bump.append(int(rng.integers(0, max(0, min(2, gaps[-1] - 2 - bump[-1])) + 1)))   # 마지막 행: 아래로
    nbk = [int(rng.integers(2, 6)) for _ in range(n)]
    levels = [yb + GRID * sum(gaps[:i]) for i in range(n)]
    pts = [(PORT_X, yb), (xl, yb)]
    for i, y in enumerate(levels):
        start_x = pts[-1][0]
        if i == n - 1:
            end_x = xl                                          # 마지막 행은 그대로 포트로
        elif i % 2 == 0:
            end_x = xr - ret[i] * GRID
        else:
            end_x = xl + ret[i] * GRID
        sign = -1 if i == n - 1 else +1
        sweep = jag_h((start_x, y), (end_x, y), sign, rng, bump[i], nbk[i])
        pts = G.join_routes(pts, sweep)
        if i + 1 < n:
            G.append_point(pts, (end_x, levels[i + 1]))
    G.append_point(pts, (PORT_X, yt))
    return G.compact_route(pts), dict(n=n, gaps=gaps, ret=ret, bump=bump, levels=[int(v) for v in levels])


def build_route(family: str, frame: G.Frame, rng, fixed=None):
    if family == "serpentine":
        return build_serp(frame, rng, fixed)
    return build_loop(frame, rng, family)


# ----------------------------------------------------------------- 게이트


def dial_str(d: dict) -> str:
    return ";".join(f"{k}={'/'.join(map(str, v)) if isinstance(v, list) else v}" for k, v in d.items())


# ----------------------------------------------------------------- 한 장
def make(index: int, seed: int, max_attempts: int = 200):
    a = assignment(index)
    r0 = np.random.default_rng(np.random.SeedSequence([seed, index, 999_999]))
    R = None
    if a["group"] != "similar":
        lo, hi = GROUPS[a["group"]]
        R = float(np.exp(r0.uniform(np.log(lo), np.log(hi))))
    why = {}
    for attempt in range(max_attempts):
        rng = np.random.default_rng(np.random.SeedSequence([seed, index, attempt]))
        try:
            if a["group"] == "similar":
                pf = p_frame_in_class(a["size_class"], rng)
                sf = similar_partner(pf, a["rel"], rng)
            else:
                fr = asym_frames(R, *GROUPS[a["group"]], rng)
                pf, sf = fr if fr else (None, None)
            if sf is None:
                why["frame"] = why.get("frame", 0) + 1; continue
            pri, dp = build_route(a["family"], pf, rng)
            fixed = dp if (a["family"] == "serpentine" and a["rel"] == "aligned") else None
            sec_left, ds = build_route(a["family"], sf, rng, fixed)
            sec = G.mirror_for_right_ports(sec_left, SIZE)
            G.validate_route(pri, SIZE); G.validate_route(sec, SIZE)
            s = G.Sample(tag=f"mix_{index:05d}", index=index, seed=seed, attempt=attempt, family=a["family"],
                         mode=f"{a['group']}/{a['rel']}", grid_um=GRID, pri_width_um=5.0, sec_width_um=5.0,
                         pri_frame=pf, sec_frame=sf, pri_route=pri, sec_route=sec, size_um=SIZE)
            lay = G.sample_layout(s)
        except (RuntimeError, ValueError):
            why["route"] = why.get("route", 0) + 1; continue
        # 실제 금속 크기(셀 기준, 검증기와 같은 정의): 각 권선이 제 벽에서 뻗은 길이(reach)와 y 범위.
        # frame이 아니라 이걸로도 0.8을 건다.
        cin = lay.cells[NET_IN][M9] | lay.cells[NET_IN][M8]; cout = lay.cells[NET_OUT][M9] | lay.cells[NET_OUT][M8]
        reach_p, reach_s = (max(x for x, _ in cin) + 1) * GRID, (N - min(x for x, _ in cout)) * GRID
        yext_p = (max(y for _, y in cin) - min(y for _, y in cin) + 1) * GRID
        yext_s = (max(y for _, y in cout) - min(y for _, y in cout) + 1) * GRID
        if a["group"] == "similar" and not (dims_ok(reach_p, reach_s) and dims_ok(yext_p, yext_s)):
            why["extent"] = why.get("extent", 0) + 1; continue
        for name, ok in (("contract", G.layout_contract_ok(lay)), ("pin_edge", pin_edge_ok(lay)),
                         ("wall", wall_ok(lay)), ("no_hole", no_hole(lay))):
            if not ok:
                why[name] = why.get(name, 0) + 1; break
        else:
            wp, hp = pf.x_right - pf.x_left, pf.y_top - pf.y_bottom
            ws, hs = sf.x_right - sf.x_left, sf.y_top - sf.y_bottom
            coincide = ""
            if a["family"] == "serpentine":                       # S 행 중 P 행과 같은 높이에 놓인 비율
                coincide = round(sum(v in set(dp["levels"]) for v in ds["levels"]) / len(ds["levels"]), 3)
            row = dict(tag=s.tag, index=index, seed=seed, attempt=attempt, **a,
                       sec_size_class=size_class_of(ws, hs),
                       target_ratio="" if R is None else round(R, 4),
                       frame_ratio=round(wp * hp / (ws * hs), 4),
                       width_ratio=round(wp / ws, 4), height_ratio=round(hp / hs, 4),
                       pri_w=wp, pri_h=hp, sec_w=ws, sec_h=hs,
                       pri_xl=pf.x_left, pri_xr=pf.x_right,
                       sec_xl_before_mirror=sf.x_left, sec_xr_before_mirror=sf.x_right,
                       sec_xl=SIZE - sf.x_right, sec_xr=SIZE - sf.x_left,
                       dx_center=(SIZE - (sf.x_left + sf.x_right) / 2) - (pf.x_left + pf.x_right) / 2,
                       pitch_in_um=hp, pitch_out_um=hs,
                       reach_pri_um=reach_p, reach_sec_um=reach_s, reach_ratio=round(reach_p / reach_s, 4),
                       yext_pri_um=yext_p, yext_sec_um=yext_s, yext_ratio=round(yext_p / yext_s, 4),
                       route_len_pri=G.route_length(pri), route_len_sec=G.route_length(sec),
                       sweeps_pri=G.sweep_count(pri), sweeps_sec=G.sweep_count(sec),
                       row_coincide_frac=coincide,
                       dials_pri=dial_str(dp), dials_sec=dial_str(ds))
            return s, row
    why.update(assignment=a, target_ratio=R)
    return None, why


def work(args):
    index, seed, outdir = args
    s, info = make(index, seed)
    if s is None:
        return index, None, info
    G.write_gds(s, Path(outdir) / f"{s.tag}.gds")                  # 셀 이름 = 파일 stem
    return index, info, None


def md5_of(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--pilot", action="store_true", help="구성별 대표 index만")
    a = ap.parse_args()
    gds_dir = Path(a.outdir) / "gds"
    if gds_dir.exists() and any(gds_dir.iterdir()):
        sys.exit(f"{gds_dir}가 비어 있지 않다. 덮어쓰지 않는다.")
    if (Path(a.outdir) / "manifest.csv").exists():
        sys.exit(f"{a.outdir}/manifest.csv가 이미 있다. 덮어쓰지 않는다.")
    gds_dir.mkdir(parents=True, exist_ok=True)
    indices = pilot_indices() if a.pilot else list(range(N_SIMILAR + 2 * N_ASYM))

    rows, fails = [], {}
    with Pool(a.workers) as pool:
        for i, info, why in pool.imap_unordered(work, [(i, a.seed, str(gds_dir)) for i in indices], chunksize=8):
            if info is None:
                fails[i] = why
            else:
                rows.append(info)
    rows.sort(key=lambda r: r["index"])
    with open(Path(a.outdir) / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    cells = Counter((r["group"], r["family"], r["size_class"], r["rel"]) for r in rows)
    json.dump(dict(n=len(rows), seed=a.seed, indices=len(indices), failed=len(fails),
                   cells={"/".join(k): v for k, v in sorted(cells.items())},
                   size_edges_um2=SIZE_EDGES.tolist(), groups=GROUPS, size_um=SIZE, grid_um=GRID,
                   generator_parts=str(GEN.relative_to(ROOT) / "generate_one_stroke_gds.py"),
                   md5=dict(gen_mix=md5_of(Path(__file__)), generate_one_stroke_gds=md5_of(GEN / "generate_one_stroke_gds.py")),
                   fail_examples={str(k): v for k, v in list(fails.items())[:20]}),
              open(Path(a.outdir) / "gen_meta.json", "w"), indent=1, default=str)
    print(f"{gds_dir}  GDS {len(rows)}장  (index {len(indices)}, 실패 {len(fails)})")
    if not a.pilot:
        assert len(rows) == N_SIMILAR + 2 * N_ASYM, f"실패한 index가 있다: {sorted(fails)[:10]}"
        assert all(v in (111, 112) for k, v in cells.items() if k[0] == "similar"), cells
        assert all(v in (333, 334) for k, v in cells.items() if k[0] != "similar"), cells


if __name__ == "__main__":
    main()
