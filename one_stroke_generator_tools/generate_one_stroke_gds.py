#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""큰 열린 루프 형태의 2층 한붓그리기 레이아웃을 GDS로 생성한다.

PRI(IN)는 M9, SEC(OUT)는 M8에 각각 하나의 끊기지 않은 선으로 배치한다.
생성된 파일은 EM 시뮬레이션 결과가 아니라 입력용 레이아웃이다.

사용 예
-------
python generate_one_stroke_gds.py --n 1000 --seed 42
python generate_one_stroke_gds.py --n 1000 --seed 42 --no-png
python generate_one_stroke_gds.py --n 100 --seed 7 --outdir ./my_dataset
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

# 홈 디렉토리가 읽기 전용인 실행 환경에서도 matplotlib 캐시 경고 없이 동작한다.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-rfic-one-stroke")
import matplotlib.pyplot as plt
import numpy as np

# GDS 레이어, 포트, VIA 기록 함수
from em_gds_contract import (  # noqa: E402
    Layout,
    M8,
    M9,
    NET_IN,
    NET_OUT,
    connectivity_ok,
    drc_ok,
    isolation_ok,
    via_ok,
    write_gds as write_em_gds,
)


# -----------------------------------------------------------------------------
# gen_dataset.py와 동일한 기본 GDS 규칙
# 좌표는 um 단위이며, GDS database unit은 0.001 um이다.
# -----------------------------------------------------------------------------
SIZE_UM = 300.0
DBU_UM = 0.001

M9_LAYER, M9_DATATYPE = 39, 60
M8_LAYER, M8_DATATYPE = 38, 40
M9_PIN_LAYER, M9_PIN_DATATYPE = 139, 0
VIA_LAYER, VIA_DATATYPE = 58, 40

PRI_COLOR = "#d12b2f"
SEC_COLOR = "#1768c5"
PORT_COLOR = "#f4f23b"

Point = tuple[float, float]


@dataclass(frozen=True)
class Frame:
    """열린 루프의 네 외곽 좌표."""

    x_left: float
    x_right: float
    y_bottom: float
    y_top: float


@dataclass
class Sample:
    """한 GDS 샘플의 형상과 기록할 메타데이터."""

    tag: str
    index: int
    seed: int
    attempt: int
    family: str
    mode: str
    grid_um: float
    pri_width_um: float
    sec_width_um: float
    pri_frame: Frame
    sec_frame: Frame
    pri_route: list[Point]
    sec_route: list[Point]


def snap(value: float, grid: float) -> float:
    """좌표를 제조 격자에 맞춘다."""

    return round(value / grid) * grid


def append_point(points: list[Point], point: Point) -> None:
    """연속 중복점을 제외하고 점 하나를 추가한다."""

    point = (float(point[0]), float(point[1]))
    if not points or point != points[-1]:
        points.append(point)


def compact_route(points: Iterable[Point]) -> list[Point]:
    """같은 직선 위의 불필요한 중간점을 제거한다."""

    compact: list[Point] = []
    for point in points:
        append_point(compact, point)

    changed = True
    while changed and len(compact) >= 3:
        changed = False
        result = [compact[0]]
        for i in range(1, len(compact) - 1):
            x0, y0 = result[-1]
            x1, y1 = compact[i]
            x2, y2 = compact[i + 1]
            if (x0 == x1 == x2) or (y0 == y1 == y2):
                changed = True
                continue
            result.append((x1, y1))
        result.append(compact[-1])
        compact = result
    return compact


def random_breaks(
    start: float,
    end: float,
    grid: float,
    count: int,
    rng: np.random.Generator,
) -> list[float]:
    """두 좌표 사이에서 격자에 맞는 서로 다른 굴곡 위치를 고른다."""

    lo, hi = sorted((start, end))
    candidates = np.arange(lo + 2 * grid, hi - 2 * grid + 0.1 * grid, grid)
    if len(candidates) == 0 or count <= 0:
        return []
    count = min(count, len(candidates))
    chosen = rng.choice(candidates, size=count, replace=False)
    return sorted(float(value) for value in chosen)


def jagged_horizontal(
    start: Point,
    end: Point,
    inward_sign: int,
    grid: float,
    rng: np.random.Generator,
    max_offset_steps: int = 2,
) -> list[Point]:
    """x 방향은 되돌아가지 않으면서 작은 직교 굴곡을 넣는다."""

    x0, y0 = start
    x1, y1 = end
    if y0 != y1:
        raise ValueError("수평 구간의 시작과 끝 y가 다릅니다.")

    forward = x1 > x0
    breaks = random_breaks(x0, x1, grid, int(rng.integers(2, 6)), rng)
    if not forward:
        breaks.reverse()

    points = [start]
    current_y = y0
    for x in breaks:
        append_point(points, (x, current_y))
        offset = int(rng.integers(0, max_offset_steps + 1)) * grid
        current_y = y0 + inward_sign * offset
        append_point(points, (x, current_y))
    append_point(points, (x1, current_y))
    append_point(points, end)
    return compact_route(points)


def jagged_vertical(
    start: Point,
    end: Point,
    inward_sign: int,
    grid: float,
    rng: np.random.Generator,
) -> list[Point]:
    """y 방향은 되돌아가지 않으면서 작은 직교 굴곡을 넣는다."""

    x0, y0 = start
    x1, y1 = end
    if x0 != x1 or y1 <= y0:
        raise ValueError("수직 구간은 같은 x에서 아래에서 위로 진행해야 합니다.")

    breaks = random_breaks(y0, y1, grid, int(rng.integers(2, 6)), rng)
    points = [start]
    current_x = x0
    for y in breaks:
        append_point(points, (current_x, y))
        offset = int(rng.integers(0, 3)) * grid
        current_x = x0 + inward_sign * offset
        append_point(points, (current_x, y))
    append_point(points, (current_x, y1))
    append_point(points, end)
    return compact_route(points)


def horizontal_in_band(
    start: Point,
    end: Point,
    y_min: float,
    y_max: float,
    grid: float,
    rng: np.random.Generator,
) -> list[Point]:
    """x는 한 방향으로만 진행하고, 지정된 y 띠 안에서는 크게 움직인다.

    기존 shallow 굴곡과 달리 루프 내부의 넓은 범위를 사용한다. x가 절대로
    되돌아가지 않기 때문에 굴곡이 깊어도 이 구간 자체는 교차하지 않는다.
    """

    x0, y0 = start
    x1, y1 = end
    if y0 != y1:
        raise ValueError("수평 자유형 구간의 시작과 끝 y가 다릅니다.")

    breaks = random_breaks(x0, x1, grid, int(rng.integers(3, 9)), rng)
    if x1 < x0:
        breaks.reverse()
    levels = np.arange(snap(y_min, grid), snap(y_max, grid) + 0.1 * grid, grid)
    if len(levels) == 0:
        levels = np.asarray([y0])

    points = [start]
    current_y = y0
    for x in breaks:
        append_point(points, (x, current_y))
        current_y = float(rng.choice(levels))
        append_point(points, (x, current_y))
    append_point(points, (x1, current_y))
    append_point(points, end)
    return compact_route(points)


def vertical_in_band(
    start: Point,
    end: Point,
    x_min: float,
    x_max: float,
    grid: float,
    rng: np.random.Generator,
) -> list[Point]:
    """y는 아래에서 위로만 진행하고, 지정된 x 띠 안에서는 크게 움직인다."""

    x0, y0 = start
    x1, y1 = end
    if x0 != x1 or y1 <= y0:
        raise ValueError("수직 자유형 구간은 같은 x에서 아래에서 위로 진행해야 합니다.")

    breaks = random_breaks(y0, y1, grid, int(rng.integers(3, 9)), rng)
    levels = np.arange(snap(x_min, grid), snap(x_max, grid) + 0.1 * grid, grid)
    if len(levels) == 0:
        levels = np.asarray([x0])

    points = [start]
    current_x = x0
    for y in breaks:
        append_point(points, (current_x, y))
        current_x = float(rng.choice(levels))
        append_point(points, (current_x, y))
    append_point(points, (current_x, y1))
    append_point(points, end)
    return compact_route(points)


def join_routes(*parts: list[Point]) -> list[Point]:
    """끝점이 같은 여러 직교 경로를 하나로 합친다."""

    joined: list[Point] = []
    for part in parts:
        for point in part:
            append_point(joined, point)
    return compact_route(joined)


def build_left_open_loop(
    frame: Frame,
    grid: float,
    rng: np.random.Generator,
) -> list[Point]:
    """왼쪽의 두 포트를 잇는 큰 열린 루프 한 줄을 만든다.

    아래쪽은 좌->우, 오른쪽은 아래->위, 위쪽은 우->좌로만 움직인다.
    서로 다른 진행 띠가 겹치지 않으므로 자기 교차가 구조적으로 방지된다.
    """

    xl, xr = frame.x_left, frame.x_right
    yb, yt = frame.y_bottom, frame.y_top
    port_x = 5.0

    lead_in = [(port_x, yb), (xl, yb)]
    bottom = jagged_horizontal((xl, yb), (xr, yb), +1, grid, rng)

    # 아래/위 굴곡과 오른쪽 굴곡이 만나는 모서리에서 충돌하지 않도록
    # 먼저 3격자만큼 직선으로 빠져나온 뒤 중앙 구간에만 굴곡을 준다.
    right_low = (xr, yb + 3 * grid)
    right_high = (xr, yt - 3 * grid)
    corner_low = [(xr, yb), right_low]
    right = jagged_vertical(right_low, right_high, -1, grid, rng)
    corner_high = [right_high, (xr, yt)]

    top = jagged_horizontal((xr, yt), (xl, yt), -1, grid, rng)
    lead_out = [(xl, yt), (port_x, yt)]
    return join_routes(lead_in, bottom, corner_low, right, corner_high, top, lead_out)


def build_deep_open_loop(
    frame: Frame,
    grid: float,
    rng: np.random.Generator,
) -> list[Point]:
    """큰 루프를 유지하면서 세 변에 깊고 비대칭적인 굴곡을 만든다."""

    xl, xr = frame.x_left, frame.x_right
    yb, yt = frame.y_bottom, frame.y_top
    port_x = 5.0
    y_mid = (yb + yt) / 2.0

    # 아래와 위의 자유형 띠 사이를 비워 서로 교차하지 않도록 한다.
    bottom_max = snap(y_mid - 3 * grid, grid)
    top_min = snap(y_mid + 3 * grid, grid)
    bottom = horizontal_in_band((xl, yb), (xr, yb), yb, bottom_max, grid, rng)
    top = horizontal_in_band((xr, yt), (xl, yt), top_min, yt, grid, rng)

    # 오른쪽 굴곡은 위/아래 자유형 띠를 지난 뒤 중앙 구간에서만 움직인다.
    right_low = (xr, bottom_max + grid)
    right_high = (xr, top_min - grid)
    right_min_x = snap(xl + 0.55 * (xr - xl), grid)
    right = vertical_in_band(right_low, right_high, right_min_x, xr, grid, rng)

    return join_routes(
        [(port_x, yb), (xl, yb)],
        bottom,
        [(xr, yb), right_low],
        right,
        [right_high, (xr, yt)],
        top,
        [(xl, yt), (port_x, yt)],
    )


def build_serpentine(
    frame: Frame,
    grid: float,
    rng: np.random.Generator,
) -> list[Point]:
    """칩 폭을 4회 또는 6회 왕복하는 다중-lobe 한붓그리기를 만든다.

    각 수평선은 서로 다른 y 띠에 있고, 연결선은 좌우 끝에만 있으므로
    별도의 탐색이나 재시도 없이도 자기 교차가 생기지 않는다.
    """

    xl, xr = frame.x_left, frame.x_right
    yb, yt = frame.y_bottom, frame.y_top
    port_x = 5.0
    sweep_count = int(rng.choice([4, 6]))

    raw_levels = np.linspace(yb, yt, sweep_count)
    y_levels = [snap(float(y), grid) for y in raw_levels]
    min_row_gap = min(b - a for a, b in zip(y_levels, y_levels[1:]))
    # 굴곡이 바로 다음 행에 닿지 않도록 항상 한 격자 이상의 간격을 남긴다.
    max_offset_steps = max(0, min(1, int(min_row_gap / grid) - 1))
    points: list[Point] = [(port_x, y_levels[0]), (xl, y_levels[0])]

    for row, y in enumerate(y_levels):
        going_right = row % 2 == 0
        start_x, end_x = (xl, xr) if going_right else (xr, xl)
        if points[-1] != (start_x, y):
            append_point(points, (start_x, y))

        # 모든 행의 작은 굴곡을 같은 방향으로 두면 이웃한 두 행이 서로를
        # 향해 움직이지 않는다. 따라서 최소 행 간격이 그대로 보존된다.
        inward_sign = +1
        sweep = jagged_horizontal(
            (start_x, y),
            (end_x, y),
            inward_sign,
            grid,
            rng,
            max_offset_steps=max_offset_steps,
        )
        points = join_routes(points, sweep)

        if row + 1 < sweep_count:
            append_point(points, (end_x, y_levels[row + 1]))

    append_point(points, (port_x, y_levels[-1]))
    return compact_route(points)


def build_family_route(
    family: str,
    frame: Frame,
    grid: float,
    rng: np.random.Generator,
) -> list[Point]:
    """선택된 형상 계열의 왼쪽 포트용 한붓그리기를 만든다."""

    if family == "large_rect":
        return build_left_open_loop(frame, grid, rng)
    if family == "deep_loop":
        return build_deep_open_loop(frame, grid, rng)
    if family == "serpentine":
        return build_serpentine(frame, grid, rng)
    raise ValueError(f"알 수 없는 형상 계열: {family}")


def mirror_for_right_ports(route: list[Point]) -> list[Point]:
    """왼쪽 포트용 루프를 x축 방향으로 뒤집어 오른쪽 포트용으로 만든다."""

    return [(SIZE_UM - x, y) for x, y in route]


def segment_intersects(a: Point, b: Point, c: Point, d: Point) -> bool:
    """두 수평/수직 폐구간이 접촉하거나 겹치는지 검사한다."""

    ax, ay = a
    bx, by = b
    cx, cy = c
    dx, dy = d
    ab_horizontal = ay == by
    cd_horizontal = cy == dy

    if ab_horizontal and cd_horizontal:
        if ay != cy:
            return False
        return max(min(ax, bx), min(cx, dx)) <= min(max(ax, bx), max(cx, dx))
    if not ab_horizontal and not cd_horizontal:
        if ax != cx:
            return False
        return max(min(ay, by), min(cy, dy)) <= min(max(ay, by), max(cy, dy))

    if ab_horizontal:
        return min(ax, bx) <= cx <= max(ax, bx) and min(cy, dy) <= ay <= max(cy, dy)
    return min(cx, dx) <= ax <= max(cx, dx) and min(ay, by) <= cy <= max(ay, by)


def validate_route(route: list[Point]) -> None:
    """직교성, 칩 경계, 자기 교차 여부를 검사한다."""

    if len(route) < 4:
        raise ValueError("경로가 너무 짧습니다.")

    for x, y in route:
        if not (0.0 <= x <= SIZE_UM and 0.0 <= y <= SIZE_UM):
            raise ValueError(f"칩 경계를 벗어난 점: {(x, y)}")

    segments = list(zip(route, route[1:]))
    for a, b in segments:
        if a == b or (a[0] != b[0] and a[1] != b[1]):
            raise ValueError(f"직교하지 않거나 길이가 0인 구간: {a} -> {b}")

    # 이웃한 두 구간은 공통 꼭짓점에서 만나는 것이 정상이다.
    for i, (a, b) in enumerate(segments):
        for j in range(i + 2, len(segments)):
            c, d = segments[j]
            if segment_intersects(a, b, c, d):
                raise ValueError(f"자기 교차: {a}->{b}, {c}->{d}")


def route_length(route: list[Point]) -> float:
    """직교 중심선의 전체 길이 [um]."""

    return sum(abs(x1 - x0) + abs(y1 - y0) for (x0, y0), (x1, y1) in zip(route, route[1:]))


def sample_frame(grid: float, rng: np.random.Generator) -> Frame:
    """300 um 영역 대부분을 쓰고 포트 중심은 y=150에 대칭인 frame."""

    # 정상 batch의 포트 pitch 60--180 um와 정확히 동일한 범위:
    # y_bottom = 150 - pitch / 2 -> 60--120 um.
    y_bottom = snap(float(rng.uniform(60.0, 121.0)), grid)
    return Frame(
        x_left=snap(float(rng.uniform(25.0, 61.0)), grid),
        x_right=snap(float(rng.uniform(239.0, 276.0)), grid),
        y_bottom=y_bottom,
        y_top=SIZE_UM - y_bottom,
    )


def shifted_frame(frame: Frame, grid: float, rng: np.random.Generator) -> Frame:
    """기준 루프를 1~2격자 이동해 부분 정렬된 두 번째 루프를 만든다."""

    dx = int(rng.choice([-2, -1, 1, 2])) * grid
    dy = int(rng.choice([-2, -1, 1, 2])) * grid
    y_bottom = snap(min(120.0, max(60.0, frame.y_bottom + dy)), grid)
    return Frame(
        x_left=snap(min(70.0, max(20.0, frame.x_left + dx)), grid),
        x_right=snap(min(280.0, max(230.0, frame.x_right + dx)), grid),
        y_bottom=y_bottom,
        y_top=SIZE_UM - y_bottom,
    )


def mirror_frame(frame: Frame) -> Frame:
    """실제 좌우 반전 뒤 원하는 위치가 되도록 frame 좌표를 변환한다."""

    return Frame(
        x_left=SIZE_UM - frame.x_right,
        x_right=SIZE_UM - frame.x_left,
        y_bottom=frame.y_bottom,
        y_top=frame.y_top,
    )


def make_sample(index: int, seed: int, attempt: int = 0) -> Sample:
    """index와 seed만으로 항상 같은 샘플을 재현한다."""

    rng = np.random.default_rng(np.random.SeedSequence([seed, index, attempt]))
    # 정상 pilot500과 동일하게 5 um 제조 격자를 고정한다.
    grid = 5.0
    pri_frame = sample_frame(grid, rng)

    # 9개마다 형상 3종 x 상대배치 3종의 모든 조합이 정확히 한 번씩 나온다.
    families = ("large_rect", "deep_loop", "serpentine")
    modes = ("aligned", "offset", "independent")
    family = families[index % len(families)]
    mode = modes[(index // len(families)) % len(modes)]

    # SEC 경로는 마지막에 좌우 반전된다. 따라서 aligned/offset은 먼저 실제
    # 목표 frame을 정한 다음 역반전하여 생성해야 물리 좌표에서 의도대로 된다.
    if mode == "aligned":
        sec_frame = mirror_frame(pri_frame)
    elif mode == "offset":
        sec_frame = mirror_frame(shifted_frame(pri_frame, grid, rng))
    else:
        sec_frame = sample_frame(grid, rng)

    pri_route = build_family_route(family, pri_frame, grid, rng)
    sec_left_route = build_family_route(family, sec_frame, grid, rng)
    sec_route = mirror_for_right_ports(sec_left_route)
    validate_route(pri_route)
    validate_route(sec_route)

    return Sample(
        tag=f"difftx_{index:05d}",
        index=index,
        seed=seed,
        attempt=attempt,
        family=family,
        mode=mode,
        grid_um=grid,
        pri_width_um=5.0,
        sec_width_um=5.0,
        pri_frame=pri_frame,
        sec_frame=sec_frame,
        pri_route=pri_route,
        sec_route=sec_route,
    )


def route_cells(route: list[Point]) -> set[tuple[int, int]]:
    """5 um 중심선 경로를 연속된 1-cell 금속 경로로 변환한다."""

    points = [(int(round(x / 5.0)), int(round(y / 5.0))) for x, y in route]
    cells: set[tuple[int, int]] = set()
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 == y1:
            cells.update((x, y0) for x in range(min(x0, x1), max(x0, x1) + 1))
        elif x0 == x1:
            cells.update((x0, y) for y in range(min(y0, y1), max(y0, y1) + 1))
        else:
            raise ValueError("직교하지 않은 route는 cell로 변환할 수 없습니다.")
    if not cells or any(not (0 <= x < 60 and 0 <= y < 60) for x, y in cells):
        raise ValueError("route cell이 60x60 canvas를 벗어났습니다.")
    return cells


def two_by_two_pin(side: str, y_center_um: float) -> list[tuple[int, int]]:
    """정상 pilot와 동일한 10x10 um 포트 pin cell 4개를 만든다."""

    center_row = int(round(y_center_um / 5.0))
    rows = (center_row - 1, center_row)
    cols = (0, 1) if side == "left" else (58, 59)
    cells = [(x, y) for x in cols for y in rows]
    if any(not (0 <= x < 60 and 0 <= y < 60) for x, y in cells):
        raise ValueError("port pin이 canvas를 벗어났습니다.")
    return cells


def sample_layout(sample: Sample) -> Layout:
    """한붓그리기 경로를 pilot500과 동일한 EM Layout 계약으로 조립한다.

    PRI는 M9에 있고, SEC body는 M8에 있다. OUT pin은 M9/M8 양쪽에 같은
    2x2 cell landing을 두고 모든 cell에 VIA8을 깔아 M8 body를 M9 port로
    확실하게 끌어올린다. 네 pin marker/label은 writer가 모두 139/0에 쓴다.
    """

    ports = {
        "IN_P": two_by_two_pin("left", sample.pri_route[0][1]),
        "IN_N": two_by_two_pin("left", sample.pri_route[-1][1]),
        "OUT_P": two_by_two_pin("right", sample.sec_route[0][1]),
        "OUT_N": two_by_two_pin("right", sample.sec_route[-1][1]),
    }
    in_pins = set(ports["IN_P"]) | set(ports["IN_N"])
    out_pins = set(ports["OUT_P"]) | set(ports["OUT_N"])
    pri = route_cells(sample.pri_route) | in_pins
    sec = route_cells(sample.sec_route) | out_pins

    layout = Layout(
        cells={
            NET_IN: {M9: pri, M8: set()},
            NET_OUT: {M9: set(out_pins), M8: sec},
        },
        vias=set(out_pins),
        ports=ports,
        floating={M9: set(), M8: set()},
        meta={"family": sample.family, "mode": sample.mode},
    )
    return layout


def layout_contract_ok(layout: Layout) -> bool:
    """정상 batch와 동일한 connectivity/isolation/via/DRC gate."""

    return (
        connectivity_ok(layout)
        and isolation_ok(layout)
        and via_ok(layout)
        and drc_ok(layout)
    )


def make_valid_sample(index: int, seed: int, max_attempts: int = 200) -> Sample:
    """형상 계열은 유지하며 모든 EM gate를 통과할 때까지 재생성한다."""

    last_error: Exception | None = None
    for attempt in range(max_attempts):
        try:
            sample = make_sample(index, seed, attempt=attempt)
            if layout_contract_ok(sample_layout(sample)):
                return sample
        except (RuntimeError, ValueError) as error:
            last_error = error
    suffix = f" 마지막 오류: {last_error}" if last_error else ""
    raise RuntimeError(f"{index}번 sample이 {max_attempts}회 안에 gate를 통과하지 못했습니다.{suffix}")


def write_gds(sample: Sample, path: Path) -> None:
    """검증된 pilot writer로 EM-equivalent GDS를 저장한다."""

    layout = sample_layout(sample)
    if not layout_contract_ok(layout):
        raise ValueError(f"{sample.tag}: EM layout gate 실패")
    write_em_gds(layout, str(path), sample.tag)


def draw_preview(sample: Sample, path: Path) -> None:
    """기존 데이터와 같은 M9/M8 2패널 PNG 미리보기를 저장한다."""

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.4), dpi=100)
    fig.suptitle(
        f"{sample.tag} | {sample.family} | {sample.mode}",
        fontsize=14,
        y=0.98,
    )

    for ax, title in zip(axes, ["M9 (PRI / IN)", "M8 (SEC / OUT)"]):
        ax.set_title(title, fontsize=13)
        ax.set_xlim(0, SIZE_UM)
        ax.set_ylim(0, SIZE_UM)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks(np.arange(0, SIZE_UM + 1, 50))
        ax.set_yticks(np.arange(0, SIZE_UM + 1, 50))

    axes[0].plot(
        *zip(*sample.pri_route),
        color=PRI_COLOR,
        linewidth=sample.pri_width_um * 1.35,
        solid_capstyle="projecting",
        solid_joinstyle="miter",
    )
    axes[1].plot(
        *zip(*sample.sec_route),
        color=SEC_COLOR,
        linewidth=sample.sec_width_um * 1.35,
        solid_capstyle="projecting",
        solid_joinstyle="miter",
    )

    ports = [
        ("IN_P", sample.pri_route[0], "left"),
        ("IN_N", sample.pri_route[-1], "left"),
        ("OUT_P", sample.sec_route[0], "right"),
        ("OUT_N", sample.sec_route[-1], "right"),
    ]
    for ax in axes:
        for name, (_, y), side in ports:
            x = 2 if side == "left" else 298
            ax.text(
                x,
                y,
                name,
                ha="left" if side == "left" else "right",
                va="center",
                fontsize=9,
                fontweight="bold",
                bbox={"facecolor": PORT_COLOR, "edgecolor": "black", "pad": 0.7},
                zorder=20,
            )

    fig.subplots_adjust(left=0.055, right=0.99, bottom=0.06, top=0.91, wspace=0.12)
    fig.savefig(path, dpi=100)
    plt.close(fig)


def manifest_row(sample: Sample) -> dict[str, object]:
    """재현과 분포 확인에 필요한 숫자를 한 행으로 만든다."""

    return {
        "tag": sample.tag,
        "index": sample.index,
        "seed": sample.seed,
        "attempt": sample.attempt,
        "family": sample.family,
        "mode": sample.mode,
        "grid_um": sample.grid_um,
        "pri_width_um": sample.pri_width_um,
        "sec_width_um": sample.sec_width_um,
        "pri_length_um": route_length(sample.pri_route),
        "sec_length_um": route_length(sample.sec_route),
        "pri_x_left": sample.pri_frame.x_left,
        "pri_x_right": sample.pri_frame.x_right,
        "pri_y_bottom": sample.pri_frame.y_bottom,
        "pri_y_top": sample.pri_frame.y_top,
        "sec_x_left_before_mirror": sample.sec_frame.x_left,
        "sec_x_right_before_mirror": sample.sec_frame.x_right,
        "sec_y_bottom": sample.sec_frame.y_bottom,
        "sec_y_top": sample.sec_frame.y_top,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="큰 열린 루프형 한붓그리기 PRI/SEC 레이아웃을 GDS로 생성",
    )
    parser.add_argument("--n", type=int, required=True, help="생성할 전체 데이터 수")
    parser.add_argument("--seed", type=int, required=True, help="전체 생성 난수 seed (0 이상의 정수)")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="출력 디렉토리. 생략하면 이 파일 옆 one_stroke_gds_seed_<seed>",
    )
    parser.add_argument("--no-png", action="store_true", help="PNG 미리보기 생성을 생략")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="같은 이름의 기존 출력 파일을 덮어씀",
    )
    args = parser.parse_args()
    if args.n <= 0:
        parser.error("--n은 1 이상이어야 합니다.")
    if args.seed < 0:
        parser.error("--seed는 0 이상이어야 합니다.")
    return args


def main() -> None:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    outdir = args.outdir or script_dir / f"one_stroke_gds_seed_{args.seed}"
    outdir.mkdir(parents=True, exist_ok=True)

    planned = [outdir / f"difftx_{i:05d}.gds" for i in range(args.n)]
    if not args.no_png:
        planned += [outdir / f"difftx_{i:05d}.png" for i in range(args.n)]
    planned += [outdir / "manifest.csv", outdir / "dataset_meta.json"]
    collisions = [path for path in planned if path.exists()]
    if collisions and not args.overwrite:
        names = ", ".join(path.name for path in collisions[:5])
        more = " ..." if len(collisions) > 5 else ""
        raise FileExistsError(
            f"기존 출력과 충돌합니다: {names}{more}\n"
            "다른 --outdir를 쓰거나, 의도한 덮어쓰기라면 --overwrite를 추가하세요."
        )

    rows: list[dict[str, object]] = []
    for index in range(args.n):
        sample = make_valid_sample(index, args.seed)
        write_gds(sample, outdir / f"{sample.tag}.gds")
        if not args.no_png:
            draw_preview(sample, outdir / f"{sample.tag}.png")
        rows.append(manifest_row(sample))

        if (index + 1) % 100 == 0 or index + 1 == args.n:
            print(f"[{index + 1:>{len(str(args.n))}}/{args.n}] 생성 완료")

    with (outdir / "manifest.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "generator": Path(__file__).name,
        "algorithm": "diverse_one_stroke_em_contract_v3",
        "count": args.n,
        "seed": args.seed,
        "size_um": SIZE_UM,
        "dbu_um": DBU_UM,
        "layers": {
            "M9_PRI": f"{M9_LAYER}/{M9_DATATYPE}",
            "M8_SEC": f"{M8_LAYER}/{M8_DATATYPE}",
            "all_pins": f"{M9_PIN_LAYER}/{M9_PIN_DATATYPE}",
            "VIA8": f"{VIA_LAYER}/{VIA_DATATYPE}",
        },
        "gds_library": "DIFF_TX.DB",
        "pixel_um": 5.0,
        "port_contract": "all ports on M9; OUT M8-M9 landing with VIA8",
        "families": ["large_rect", "deep_loop", "serpentine"],
        "modes": ["aligned", "offset", "independent"],
        "png_generated": not args.no_png,
    }
    with (outdir / "dataset_meta.json").open("w", encoding="utf-8") as file:
        json.dump(metadata, file, ensure_ascii=False, indent=2)
        file.write("\n")

    print(f"출력: {outdir.resolve()}")


if __name__ == "__main__":
    main()
