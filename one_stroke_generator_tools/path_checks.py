"""5 um 금속 셀의 단순 경로 검사. 10 um 포트 패드는 각각 한 노드로 본다."""


def simple_path_ok(cells, first_pad, last_pad):
    cells, first_pad, last_pad = set(cells), set(first_pad), set(last_pad)
    if not first_pad or not last_pad or first_pad & last_pad:
        return False
    if not (first_pad | last_pad) <= cells:
        return False

    def node(p):
        return "first" if p in first_pad else "last" if p in last_pad else p

    graph = {node(p): set() for p in cells}
    for x, y in cells:
        a = node((x, y))
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            q = (x + dx, y + dy)
            if q in cells and node(q) != a:
                graph[a].add(node(q))
        # 정상 L자 코너 외에, 대각선으로만 만나는 금속 구간도 제외한다.
        for dx, dy in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            q = (x + dx, y + dy)
            if q in cells and node(q) != a and (x + dx, y) not in cells and (x, y + dy) not in cells:
                return False
    if any(len(edges) != (1 if p in ("first", "last") else 2) for p, edges in graph.items()):
        return False
    seen, stack = set(), ["first"]
    while stack:
        p = stack.pop()
        if p not in seen:
            seen.add(p)
            stack.extend(graph[p] - seen)
    return len(seen) == len(graph)
