#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GDS 레이어, 포트, VIA 기록과 레이아웃 검증 함수를 제공한다.

GDS 기록 규칙은 pilot500의 ``pilotgen.gdsio.write_gds``와 동일하다.
형상 생성은 ``generate_one_stroke_gds.py``에서 수행한다.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

import gdstk


# ---------------------------------------------------------------------------
# GDS 출력 규격
# ---------------------------------------------------------------------------
GRID_UM = 5.0
N = 60
SIZE_UM = 300.0

L_M9 = (39, 60)
L_M8 = (38, 40)
L_VIA = (58, 40)
L_PIN = (139, 0)
LIB_NAME = "DIFF_TX.DB"

M9, M8 = 0, 1
NET_IN, NET_OUT = 1, 2
PORT_NAMES = ("IN_P", "IN_N", "OUT_P", "OUT_N")

VIA_SIZE = 0.36
VIA_PX, VIA_PY = 0.92, 0.70
VIA_ENC = 0.30


@dataclass
class Layout:
    """GDS로 기록하기 전의 5 um cell 표현."""

    # net -> layer -> {(ix, iy)}
    cells: dict[int, dict[int, set[tuple[int, int]]]]
    vias: set[tuple[int, int]]
    ports: dict[str, list[tuple[int, int]]]
    floating: dict[int, set[tuple[int, int]]] = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def all_metal(self, layer: int) -> set[tuple[int, int]]:
        metal: set[tuple[int, int]] = set()
        for net in (NET_IN, NET_OUT):
            metal |= self.cells.get(net, {}).get(layer, set())
        metal |= self.floating.get(layer, set())
        return metal


def _neighbours_3d(cell, layer, vias):
    ix, iy = cell
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        yield (ix + dx, iy + dy), layer
    if cell in vias:
        yield cell, 1 - layer


def connectivity_ok(layout: Layout) -> bool:
    """IN과 OUT이 각각 P에서 N까지 DC 연결되는지 검사한다."""

    for net, names in (
        (NET_IN, ("IN_P", "IN_N")),
        (NET_OUT, ("OUT_P", "OUT_N")),
    ):
        occupied = {
            layer: set(cells)
            for layer, cells in layout.cells.get(net, {}).items()
        }
        seeds = [
            (cell, M9)
            for cell in layout.ports[names[0]]
            if cell in occupied.get(M9, set())
        ]
        if not seeds:
            return False

        seen = set(seeds)
        stack = list(seeds)
        while stack:
            cell, layer = stack.pop()
            for next_cell, next_layer in _neighbours_3d(cell, layer, layout.vias):
                if (next_cell, next_layer) in seen:
                    continue
                if next_cell not in occupied.get(next_layer, set()):
                    continue
                seen.add((next_cell, next_layer))
                stack.append((next_cell, next_layer))

        if not any((cell, M9) in seen for cell in layout.ports[names[1]]):
            return False
    return True


def isolation_ok(layout: Layout) -> bool:
    """PRI와 SEC가 같은 금속층에서 직접 닿지 않는지 검사한다."""

    for layer in (M9, M8):
        pri = layout.cells.get(NET_IN, {}).get(layer, set())
        sec = layout.cells.get(NET_OUT, {}).get(layer, set())
        if pri & sec:
            return False
    return True


def via_ok(layout: Layout) -> bool:
    """모든 VIA가 같은 net의 M9/M8 metal 사이에 있는지 검사한다."""

    for cell in layout.vias:
        owner = None
        for net in (NET_IN, NET_OUT):
            if cell in layout.cells.get(net, {}).get(M9, set()):
                owner = net
        if owner is None:
            return False
        if cell not in layout.cells.get(owner, {}).get(M8, set()):
            return False
    return True


def _has_solid_block(cells: set[tuple[int, int]], size: int = 3) -> bool:
    for ix, iy in cells:
        if all(
            (ix + dx, iy + dy) in cells
            for dx in range(size)
            for dy in range(size)
        ):
            return True
    return False


def _different_net_spacing_ok(a: set, b: set, min_gap_cells: int = 2) -> bool:
    """서로 다른 net 사이에 최소 한 cell의 빈 공간이 있는지 검사한다."""

    if not a or not b:
        return True
    radius = min_gap_cells - 1
    halo = {
        (ix + dx, iy + dy)
        for ix, iy in a
        for dx in range(-radius, radius + 1)
        for dy in range(-radius, radius + 1)
    }
    return not (halo & b)


def drc_ok(layout: Layout) -> bool:
    """10 um 최대 폭과 서로 다른 net의 최소 간격을 검사한다."""

    for layer in (M9, M8):
        pri = layout.cells.get(NET_IN, {}).get(layer, set())
        sec = layout.cells.get(NET_OUT, {}).get(layer, set())
        floating = layout.floating.get(layer, set())
        if _has_solid_block(pri | sec | floating):
            return False
        if not _different_net_spacing_ok(pri, sec):
            return False
        if not _different_net_spacing_ok(pri, floating):
            return False
        if not _different_net_spacing_ok(sec, floating):
            return False
    return True


def merge_cells(cells: set[tuple[int, int]]) -> list[tuple[int, int, int, int]]:
    """5 um cell 집합을 겹치지 않는 직사각형 목록으로 병합한다."""

    remaining = set(cells)
    rectangles = []
    while remaining:
        ix, iy = min(remaining, key=lambda cell: (cell[1], cell[0]))
        x1 = ix
        while (x1 + 1, iy) in remaining:
            x1 += 1
        y1 = iy
        while all((x, y1 + 1) in remaining for x in range(ix, x1 + 1)):
            y1 += 1
        for y in range(iy, y1 + 1):
            for x in range(ix, x1 + 1):
                remaining.discard((x, y))
        rectangles.append((ix, iy, x1 + 1, y1 + 1))
    return rectangles


def _via_array(x0: float, y0: float, width: float, height: float):
    """하나의 metal cell 안에 들어가는 VIA8 cut의 좌하단 좌표."""

    nx = int(math.floor((width - 2 * VIA_ENC - VIA_SIZE) / VIA_PX)) + 1
    ny = int(math.floor((height - 2 * VIA_ENC - VIA_SIZE) / VIA_PY)) + 1
    if nx < 1 or ny < 1:
        return []
    margin_x = (width - ((nx - 1) * VIA_PX + VIA_SIZE)) / 2.0
    margin_y = (height - ((ny - 1) * VIA_PY + VIA_SIZE)) / 2.0
    return [
        (x0 + margin_x + i * VIA_PX, y0 + margin_y + j * VIA_PY)
        for j in range(ny)
        for i in range(nx)
    ]


def _via_cuts_for_cell(ix: int, iy: int):
    return _via_array(ix * GRID_UM, iy * GRID_UM, GRID_UM, GRID_UM)


def _label_origin(cells: list[tuple[int, int]]) -> tuple[float, float]:
    ix = min(cell[0] for cell in cells)
    iy = min(cell[1] for cell in cells)
    return (ix * GRID_UM + GRID_UM / 2.0, iy * GRID_UM + GRID_UM / 2.0)


def write_gds(layout: Layout, path: str, cell_name: str) -> None:
    """검증된 EM 계약으로 GDS 하나를 기록한다."""

    library = gdstk.Library(name=LIB_NAME, unit=1e-6, precision=1e-9)
    cell = library.new_cell(cell_name)

    for layer_index, (layer, datatype) in ((M9, L_M9), (M8, L_M8)):
        for ix0, iy0, ix1, iy1 in merge_cells(layout.all_metal(layer_index)):
            cell.add(
                gdstk.rectangle(
                    (ix0 * GRID_UM, iy0 * GRID_UM),
                    (ix1 * GRID_UM, iy1 * GRID_UM),
                    layer=layer,
                    datatype=datatype,
                )
            )

    for via_cell in sorted(layout.vias):
        for x, y in _via_cuts_for_cell(*via_cell):
            cell.add(
                gdstk.rectangle(
                    (x, y),
                    (x + VIA_SIZE, y + VIA_SIZE),
                    layer=L_VIA[0],
                    datatype=L_VIA[1],
                )
            )

    # 포트 marker와 label은 네 개 모두 M9 pin layer 139/0이다.
    for name in PORT_NAMES:
        cells = layout.ports[name]
        ix0 = min(cell_[0] for cell_ in cells)
        ix1 = max(cell_[0] for cell_ in cells) + 1
        iy0 = min(cell_[1] for cell_ in cells)
        iy1 = max(cell_[1] for cell_ in cells) + 1
        cell.add(
            gdstk.rectangle(
                (ix0 * GRID_UM, iy0 * GRID_UM),
                (ix1 * GRID_UM, iy1 * GRID_UM),
                layer=L_PIN[0],
                datatype=L_PIN[1],
            )
        )
        cell.add(
            gdstk.Label(
                name,
                _label_origin(cells),
                anchor="nw",
                magnification=1.0,
                rotation=0.0,
                layer=L_PIN[0],
                texttype=L_PIN[1],
            )
        )

    library.write_gds(path)
