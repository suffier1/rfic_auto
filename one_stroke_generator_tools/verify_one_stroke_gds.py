#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""one-stroke GDS 디렉토리의 EM 규격과 전기적 연결을 전수검사한다.

사용 예
-------
python verify_one_stroke_gds.py ./one_stroke_gds_seed_8_corrected --expected-count 1000
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from pathlib import Path
import sys

import gdstk

from em_gds_contract import (
    GRID_UM,
    LIB_NAME,
    L_M8,
    L_M9,
    L_PIN,
    L_VIA,
    N,
    PORT_NAMES,
)


ALLOWED_LAYERS = {L_M9, L_M8, L_VIA, L_PIN}


def polygon_cells(polygon) -> set[tuple[int, int]]:
    """5 um 직사각형 metal을 cell 집합으로 되돌린다."""

    (x0, y0), (x1, y1) = polygon.bounding_box()
    ix0, iy0 = int(round(x0 / GRID_UM)), int(round(y0 / GRID_UM))
    ix1, iy1 = int(round(x1 / GRID_UM)), int(round(y1 / GRID_UM))
    return {(ix, iy) for ix in range(ix0, ix1) for iy in range(iy0, iy1)}


def reachable(
    start: tuple[int, int, int],
    metal: dict[int, set[tuple[int, int]]],
    via_cells: set[tuple[int, int]],
) -> set[tuple[int, int, int]]:
    """M9/M8와 VIA를 따라 start에서 DC로 닿는 모든 cell을 찾는다."""

    seen = {start}
    queue = deque([start])
    while queue:
        layer, ix, iy = queue.popleft()
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            candidate = (layer, ix + dx, iy + dy)
            if candidate in seen:
                continue
            if (candidate[1], candidate[2]) in metal[layer]:
                seen.add(candidate)
                queue.append(candidate)
        if (ix, iy) in via_cells:
            candidate = (1 - layer, ix, iy)
            if candidate not in seen and (ix, iy) in metal[1 - layer]:
                seen.add(candidate)
                queue.append(candidate)
    return seen


def verify_one(path: Path) -> list[str]:
    """파일 하나를 검사하고 실패 사유 목록을 반환한다."""

    errors: list[str] = []
    try:
        library = gdstk.read_gds(path)
    except Exception as error:
        return [f"GDS 읽기 실패: {error}"]

    top = library.top_level()
    if library.name != LIB_NAME:
        errors.append(f"library={library.name!r}, expected={LIB_NAME!r}")
    if library.unit != 1e-6 or library.precision != 1e-9:
        errors.append(f"unit/precision={library.unit}/{library.precision}")
    if len(top) != 1:
        return errors + [f"top cell 수={len(top)}, expected=1"]

    cell = top[0]
    if cell.name != path.stem:
        errors.append(f"cell name={cell.name!r}, file stem={path.stem!r}")
    if cell.references or cell.paths:
        errors.append("flat GDS가 아님: reference 또는 PATH 존재")
    if not all(len(polygon.points) == 4 for polygon in cell.polygons):
        errors.append("4점 직사각형이 아닌 polygon 존재")

    layer_counts = Counter((p.layer, p.datatype) for p in cell.polygons)
    if not set(layer_counts) <= ALLOWED_LAYERS:
        errors.append(f"허용되지 않은 layer 존재: {set(layer_counts) - ALLOWED_LAYERS}")
    for required in (L_M9, L_M8, L_PIN, L_VIA):
        if required not in layer_counts:
            errors.append(f"필수 layer 없음: {required}")

    pins = [p for p in cell.polygons if (p.layer, p.datatype) == L_PIN]
    labels = {label.text: label for label in cell.labels}
    if len(pins) != 4:
        errors.append(f"pin marker 수={len(pins)}, expected=4")
    if set(labels) != set(PORT_NAMES):
        errors.append(f"label 이름={sorted(labels)}, expected={sorted(PORT_NAMES)}")

    for label in cell.labels:
        if (label.layer, label.texttype) != L_PIN or label.anchor != "nw":
            errors.append(
                f"{label.text}: layer/texttype/anchor="
                f"{label.layer}/{label.texttype}/{label.anchor}"
            )

    # label origin은 각 10x10 um pin의 좌하단 5 um cell 중심이다.
    port_cells: dict[str, tuple[int, int]] = {}
    for name, label in labels.items():
        ix = int((label.origin[0] - GRID_UM / 2) / GRID_UM + 0.5)
        iy = int((label.origin[1] - GRID_UM / 2) / GRID_UM + 0.5)
        port_cells[name] = (ix, iy)
        expected_x = 0 if name.startswith("IN") else 58
        if ix != expected_x:
            errors.append(f"{name}: label x cell={ix}, expected={expected_x}")

    if len(port_cells) == 4:
        pitch_in = abs(port_cells["IN_N"][1] - port_cells["IN_P"][1]) * GRID_UM
        pitch_out = abs(port_cells["OUT_N"][1] - port_cells["OUT_P"][1]) * GRID_UM
        for name, pitch in (("IN", pitch_in), ("OUT", pitch_out)):
            if pitch not in range(60, 181, 10):
                errors.append(f"{name} pitch={pitch} um, expected=60..180 step 10")
        for side in ("IN", "OUT"):
            y_sum = (
                port_cells[f"{side}_P"][1]
                + port_cells[f"{side}_N"][1]
                + 2
            ) * GRID_UM
            if y_sum != 300.0:
                errors.append(f"{side} pin 중심이 y=150 대칭이 아님")

    metal = {0: set(), 1: set()}
    via_cells: set[tuple[int, int]] = set()
    for polygon in cell.polygons:
        layer_key = (polygon.layer, polygon.datatype)
        if layer_key == L_M9:
            metal[0] |= polygon_cells(polygon)
        elif layer_key == L_M8:
            metal[1] |= polygon_cells(polygon)
        elif layer_key == L_VIA:
            (x0, y0), _ = polygon.bounding_box()
            via_cells.add((int(x0 / GRID_UM), int(y0 / GRID_UM)))

    if len(via_cells) != 8:
        errors.append(f"VIA landing cell 수={len(via_cells)}, expected=8")
    if layer_counts.get(L_VIA, 0) != 240:
        errors.append(f"VIA8 cut 수={layer_counts.get(L_VIA, 0)}, expected=240")

    if len(port_cells) == 4:
        starts = {
            name: (0, *cell_xy)
            for name, cell_xy in port_cells.items()
        }
        for name, start in starts.items():
            if (start[1], start[2]) not in metal[0]:
                errors.append(f"{name}: M9 metal에 붙어 있지 않음")

        if not errors:
            pri_seen = reachable(starts["IN_P"], metal, via_cells)
            sec_seen = reachable(starts["OUT_P"], metal, via_cells)
            if starts["IN_N"] not in pri_seen:
                errors.append("IN_P와 IN_N이 DC로 연결되지 않음")
            if starts["OUT_N"] not in sec_seen:
                errors.append("OUT_P와 OUT_N이 DC로 연결되지 않음")
            if starts["OUT_P"] in pri_seen or starts["OUT_N"] in pri_seen:
                errors.append("PRI와 SEC가 DC short됨")
            if starts["IN_P"] in sec_seen or starts["IN_N"] in sec_seen:
                errors.append("SEC와 PRI가 DC short됨")

    return errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="one-stroke GDS EM 규격 전수검사")
    parser.add_argument("directory", type=Path, help="difftx_*.gds가 있는 디렉토리")
    parser.add_argument("--expected-count", type=int, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    files = sorted(args.directory.glob("difftx_*.gds"))
    if args.expected_count is not None and len(files) != args.expected_count:
        print(f"FAIL: GDS 수={len(files)}, expected={args.expected_count}")
        return 1
    if not files:
        print("FAIL: difftx_*.gds가 없습니다.")
        return 1

    failed = []
    for index, path in enumerate(files, 1):
        errors = verify_one(path)
        if errors:
            failed.append((path, errors))
        if index % 100 == 0 or index == len(files):
            print(f"[{index}/{len(files)}] 검사 완료")

    if failed:
        print(f"\nFAIL: {len(failed)}/{len(files)} files")
        for path, errors in failed[:10]:
            print(f"- {path.name}: {'; '.join(errors)}")
        return 1

    print(f"\nPASS: {len(files)}/{len(files)} GDS 모두 EM 규격과 DC 연결 정상")
    return 0


if __name__ == "__main__":
    sys.exit(main())
