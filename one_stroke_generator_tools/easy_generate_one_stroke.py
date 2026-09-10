#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""설정값에 따라 one-stroke GDS를 생성하고 전수검사한다."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


# ===========================================================================
# 생성 설정
# ===========================================================================
TOTAL_COUNT = 1000
RANDOM_SEED = 8
OUTPUT_FOLDER_NAME = "one_stroke_gds_seed_8_new"
MAKE_PNG = False
FAMILY = "all"  # all / large_rect / deep_loop / serpentine
START_INDEX = 0  # 같은 seed로 분할 생성할 때 이전 구간 다음 번호
SERPENTINE_ENVELOPE = "bounded"  # bounded / expanded / mixed
# ===========================================================================


def main() -> None:
    parser = argparse.ArgumentParser(description="GDS 생성과 전수검사를 순서대로 실행")
    parser.add_argument("--n", type=int, default=TOTAL_COUNT, help="이 batch에서 생성할 개수")
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--family", choices=("all", "large_rect", "deep_loop", "serpentine"), default=FAMILY)
    parser.add_argument("--start-index", type=int, default=START_INDEX)
    parser.add_argument("--serpentine-envelope", choices=("bounded", "expanded", "mixed"), default=SERPENTINE_ENVELOPE)
    parser.add_argument("--outdir", type=Path, default=None, help="지정하면 현재 작업 디렉터리 기준 출력 경로")
    parser.add_argument("--png", action=argparse.BooleanOptionalAction, default=MAKE_PNG)
    args = parser.parse_args()
    script_dir = Path(__file__).resolve().parent
    generator = script_dir / "generate_one_stroke_gds.py"
    verifier = script_dir / "verify_one_stroke_gds.py"
    output_dir = (args.outdir or script_dir / OUTPUT_FOLDER_NAME).resolve()
    if args.n <= 0 or args.seed < 0 or args.start_index < 0:
        parser.error("n은 1 이상, seed와 start-index는 0 이상이어야 합니다.")
    if output_dir.exists():
        raise FileExistsError(
            "같은 출력 디렉터리가 이미 있습니다. 기존 파일을 자동 삭제하지 않습니다.\n"
            "OUTPUT_FOLDER_NAME을 새 이름으로 바꾸고 다시 실행하세요.\n"
            f"directory: {output_dir}"
        )

    generate_command = [
        sys.executable,
        str(generator),
        "--n",
        str(args.n),
        "--seed",
        str(args.seed),
        "--family",
        args.family,
        "--start-index",
        str(args.start_index),
        "--serpentine-envelope",
        args.serpentine_envelope,
        "--outdir",
        str(output_dir),
    ]
    if not args.png:
        generate_command.append("--no-png")

    print("[1/2] GDS 생성", flush=True)
    subprocess.run(generate_command, check=True)

    print("[2/2] 전체 GDS 검증", flush=True)
    subprocess.run(
        [
            sys.executable,
            str(verifier),
            str(output_dir),
            "--expected-count",
            str(args.n),
            "--report",
            str(output_dir / "verification.json"),
        ],
        check=True,
    )

    print("\n완료")
    print(f"GDS directory: {output_dir}")


if __name__ == "__main__":
    main()
