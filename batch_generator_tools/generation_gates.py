"""배치 생성 중에 쓰는 기하 게이트. 생성기(gen_asym, gen_mix)가 떨어진 경로를 버리고 다시 뽑을 때 쓴다.

- pin_edge_ok  같은 쪽 두 핀이 가장자리 열의 연속 금속으로 이어지면 False (Momentum: 한 에지에 pin 하나)
- wall_ok      좌우 벽 두 열(0~1, 58~59)에 단자가 정확히 2개씩
- not_shorted  권선이 제 금속으로 50칸 넘는 구멍을 감싸지 않는다 (gen_asym)
- no_hole      구멍이 하나도 없다, 크기 무관 (gen_mix)
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

# em_gds_contract와 같은 값. gdstk 없이도 이 파일을 쓸 수 있게 직접 둔다.
N = 60
M9, M8 = 0, 1
NET_IN, NET_OUT = 1, 2


def self_closed(cells) -> bool:
    """셀 집합이 50칸 넘는 구멍을 감싸는가. (개인 코드의 serp 판정과 같은 정의 — 그쪽은 gdstk를
    끌어와서 torch 환경에서 못 쓰므로 여기 따로 둔다)"""
    m = np.zeros((N, N), bool)
    for x, y in cells:
        if 0 <= x < N and 0 <= y < N:
            m[y, x] = True
    return int((ndimage.binary_fill_holes(m) & ~m).sum()) > 50


EDGE = {"IN": (0, ("IN_P", "IN_N")), "OUT": (N - 1, ("OUT_P", "OUT_N"))}


def pin_edge_ok(lay) -> bool:
    """같은 쪽 두 핀 사이의 가장자리 열이 금속으로 끊김 없이 이어지면 False."""
    m9 = lay.all_metal(M9)
    for side, (col, (p, n)) in EDGE.items():
        rp = max(c[1] for c in lay.ports[p])            # P 핀 윗끝 행
        rn = min(c[1] for c in lay.ports[n])            # N 핀 아랫끝 행
        lo, hi = min(rp, rn), max(rp, rn)
        if hi - lo <= 1:
            return False                                # 두 핀이 붙어 있다
        if all((col, r) in m9 for r in range(lo + 1, hi)):
            return False                                # 가장자리 열이 두 핀을 잇는다
    return True


def not_shorted(lay) -> bool:
    """어느 net도 50칸 넘는 구멍을 감싸지 않는다 (옛 기준 — gen_asym.py 전용)."""
    return not any(self_closed(lay.cells[k][M9] | lay.cells[k][M8]) for k in (NET_IN, NET_OUT))


SIDES = (("IN", 1, (0, 1)), ("OUT", -1, (N - 2, N - 1)))


def runs_of(rows: np.ndarray) -> list[list[int]]:
    out = []
    for r in np.nonzero(rows)[0].tolist():
        if out and r == out[-1][-1] + 1:
            out[-1].append(r)
        else:
            out.append([r])
    return out


def derive(m9: np.ndarray, snap: bool = False):
    """(포트 지도, 사유). 규격을 못 맞추면 (None, 사유).

    snap=True면 P/N 핀을 캔버스 중심(150 um)에 대칭인 가장 가까운 자리로 옮긴다.
    배달분 규약은 150/150 파일이 대칭인데 무조건부 모델은 28%를 어긋나게 그린다 —
    대부분 한두 칸 차이라 버리는 대신 맞춰 준다. 옮긴 자리는 m9에 금속을 채워
    권선이 새 핀에 닿게 한다(호출한 쪽이 m9를 넘겨받아 반영해야 한다).
    """
    p = np.zeros((N, N), np.int8)
    for side, val, cols in SIDES:
        rs = runs_of(m9[:, cols[0]] | m9[:, cols[1]])
        if len(rs) != 2:
            return None, f"{side} 단자 {len(rs)}개 (2개여야 함)"
        lo, hi = rs[0], rs[-1]
        if lo[-1] >= hi[0]:
            return None, f"{side} 단자가 붙어 있음"
        w = min(2, len(lo), len(hi))                     # 핀 높이 (양쪽 같게)
        if snap:
            rp = int(round((lo[0] + (N - w - hi[0])) / 2))   # 대칭 조건: rn = N - w - rp
            rp = max(0, min(rp, N // 2 - w))
            rows_p = list(range(rp, rp + w))
            rows_n = list(range(N - w - rp, N - rp))
        else:
            rows_p = lo[:2] if len(lo) > 2 else lo
            rows_n = hi[:2] if len(hi) > 2 else hi
        for rows in (rows_p, rows_n):
            for r in rows:
                for c in cols:                           # x로 정확히 2칸
                    p[r, c] = val
    return p, ""


def wall_ok(layout) -> bool:
    """좌우 벽 두 열에 단자가 정확히 2개씩 — 벽을 따라 두 핀을 잇는 금속이 없다."""
    m9 = np.zeros((60, 60), bool)
    for x, y in layout.all_metal(M9):
        m9[y, x] = True
    return derive(m9)[0] is not None


def hole_cells(lay, net) -> int:
    m = np.zeros((N, N), bool)
    for x, y in lay.cells[net][M9] | lay.cells[net][M8]:
        m[y, x] = True
    return int((ndimage.binary_fill_holes(m) & ~m).sum())


def no_hole(lay) -> bool:
    return hole_cells(lay, NET_IN) == 0 and hole_cells(lay, NET_OUT) == 0
