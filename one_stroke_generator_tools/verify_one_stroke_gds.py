#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""one-stroke GDS 디렉토리의 EM 규격과 전기적 연결을 전수검사한다.

사용 예
-------
python verify_one_stroke_gds.py ./one_stroke_gds_seed_8_corrected --expected-count 1000
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, deque
from pathlib import Path
import sys

import gdstk
import numpy as np

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


def verify_one(path: Path, size_um: int = 300, strict_path: bool = False) -> list[str]:
    """파일 하나를 검사하고 실패 사유 목록을 반환한다."""

    errors: list[str] = []
    if size_um not in (200, 300):
        return ["지원하지 않는 canvas 크기"]
    n = int(size_um / GRID_UM)
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
    for p in cell.polygons:
        lo, hi = np.asarray(p.bounding_box())
        if abs(p.area() - np.prod(hi - lo)) > 1e-7:
            errors.append("직사각형이 아닌 polygon")
        if np.any(lo < 0) or np.any(hi > size_um):
            errors.append("canvas 경계 초과")
        if (p.layer, p.datatype) != L_VIA:
            values = np.concatenate((lo, hi)) / GRID_UM
            if not np.allclose(values, np.rint(values), atol=1e-7, rtol=0):
                errors.append("5 um 격자 불일치")

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
    if set(labels) != set(PORT_NAMES) or len(cell.labels) != 4:
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
        expected_x = 0 if name.startswith("IN") else n - 2
        if ix != expected_x:
            errors.append(f"{name}: label x cell={ix}, expected={expected_x}")

    if set(port_cells) == set(PORT_NAMES):
        pitch_in = abs(port_cells["IN_N"][1] - port_cells["IN_P"][1]) * GRID_UM
        pitch_out = abs(port_cells["OUT_N"][1] - port_cells["OUT_P"][1]) * GRID_UM
        for name, pitch in (("IN", pitch_in), ("OUT", pitch_out)):
            allowed = range(60, 181, 10) if size_um == 300 else range(40, 121, 10)
            if pitch not in allowed:
                errors.append(f"{name} pitch={pitch} um, 허용 범위 불일치")
        for side in ("IN", "OUT"):
            y_sum = (
                port_cells[f"{side}_P"][1]
                + port_cells[f"{side}_N"][1]
                + 2
            ) * GRID_UM
            if y_sum != size_um:
                errors.append(f"{side} pin 중심이 y={size_um / 2:g} 대칭이 아님")

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

    if set(port_cells) == set(PORT_NAMES):
        pads = {name: {(x + dx, y + dy) for dx in (0, 1) for dy in (0, 1)}
                for name, (x, y) in port_cells.items()}
        actual_pins = [polygon_cells(p) for p in pins]
        if any(actual_pins.count(pad) != 1 for pad in pads.values()):
            errors.append("10x10 um pin marker와 label 위치 불일치")
        out_pads = pads["OUT_P"] | pads["OUT_N"]
        if via_cells != out_pads or not out_pads <= metal[0] & metal[1]:
            errors.append("OUT landing/VIA 위치 불일치")
        # 5x5 um landing cell마다 5열 x 6행 cut, 기존 공정 수치 고정.
        expected_vias = Counter()
        for x, y in out_pads:
            for i in range(5):
                for j in range(6):
                    x0, y0 = x * 5 + .48 + i * .92, y * 5 + .57 + j * .70
                    expected_vias[tuple(round(v * 1000) for v in (x0, y0, x0 + .36, y0 + .36))] += 1
        actual_vias = Counter(tuple(round(float(v) * 1000) for point in p.bounding_box() for v in point)
                              for p in cell.polygons if (p.layer, p.datatype) == L_VIA)
        if actual_vias != expected_vias:
            errors.append("VIA cut 크기/간격/배열 불일치")
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
            all_nodes = {(layer, x, y) for layer, cells in metal.items() for x, y in cells}
            if pri_seen | sec_seen != all_nodes:
                errors.append("포트와 연결되지 않은 floating metal")
            if any(layer != 0 for layer, _, _ in pri_seen):
                errors.append("PRI가 M9 밖에 존재")
            if {(x, y) for layer, x, y in sec_seen if layer == 0} != out_pads:
                errors.append("SEC의 M9 금속이 OUT landing 외부에 존재")
            pri = {(x, y) for layer, x, y in pri_seen if layer == 0}
            if any((x + dx, y + dy) in out_pads for x, y in pri
                   for dx in (-1, 0, 1) for dy in (-1, 0, 1)):
                errors.append("M9의 서로 다른 net 간격 부족")
            if strict_path:
                from path_checks import simple_path_ok
                for name, cells in (("IN", pri), ("OUT", metal[1])):
                    if not simple_path_ok(cells, pads[name + "_P"], pads[name + "_N"]):
                        errors.append(f"{name}: 분기·폐회로·비인접 접촉 또는 경로 단절")

    return errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="one-stroke GDS EM 규격 전수검사")
    parser.add_argument("directory", type=Path, help="difftx_*.gds가 있는 디렉토리")
    parser.add_argument("--expected-count", type=int, default=None)
    parser.add_argument("--size-um", type=int, choices=(200, 300), default=None, help="생략하면 dataset_meta.json에서 읽음; 구형 데이터 기본 300")
    parser.add_argument("--strict-path", action="store_true", help="포트 패드 포함 단순 경로 검사; v5에서는 자동 적용")
    parser.add_argument("--check-png", action="store_true", help="PNG 금속 영역과 GDS의 픽셀 일치 검사")
    parser.add_argument("--report", type=Path, default=None, help="검사 결과 JSON 저장 경로")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    files = sorted(args.directory.glob("difftx_*.gds"))
    meta_path = args.directory / "dataset_meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    size_um = args.size_um if args.size_um is not None else meta.get("size_um", 300)
    strict = args.strict_path or meta.get("serpentine_envelope") == "freeform"
    count_ok = bool(files) and (args.expected_count is None or len(files) == args.expected_count)

    failed = []
    for index, path in enumerate(files, 1):
        errors = verify_one(path, size_um=size_um, strict_path=strict)
        if args.check_png:
            from render_gds import preview_matches
            try:
                if not preview_matches(path, path.with_suffix(".png"), size_um):
                    errors.append("PNG/GDS 금속 영역 불일치")
            except Exception as error:
                errors.append(f"PNG 검사 실패: {error}")
        if errors:
            failed.append((path, errors))
        if index % 100 == 0 or index == len(files):
            print(f"[{index}/{len(files)}] 검사 완료")

    if args.report:
        report = {
            "passed": count_ok and not failed,
            "count": len(files),
            "expected_count": args.expected_count,
            "size_um": size_um,
            "strict_path": strict,
            "png_checked": args.check_png,
            "failed": [{"file": p.name, "errors": errors} for p, errors in failed],
            "scope": "GDS contract and DC connectivity; not foundry DRC or EM validation",
        }
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if not count_ok:
        print(f"FAIL: GDS 수={len(files)}, expected={args.expected_count}; 빈 디렉터리는 통과하지 않습니다.")
        return 1
    if failed:
        print(f"\nFAIL: {len(failed)}/{len(files)} files")
        for path, errors in failed[:10]:
            print(f"- {path.name}: {'; '.join(errors)}")
        return 1

    print(f"\nPASS: {len(files)}/{len(files)} GDS 모두 EM 규격과 DC 연결 정상")
    return 0


if __name__ == "__main__":
    sys.exit(main())
