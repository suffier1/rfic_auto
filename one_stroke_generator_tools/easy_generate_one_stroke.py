#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""설정값에 따라 one-stroke GDS를 생성하고 전수검사한다."""

from __future__ import annotations

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
# ===========================================================================


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    generator = script_dir / "generate_one_stroke_gds.py"
    verifier = script_dir / "verify_one_stroke_gds.py"
    output_dir = script_dir / OUTPUT_FOLDER_NAME
    if TOTAL_COUNT <= 0:
        raise ValueError("TOTAL_COUNT는 1 이상이어야 합니다.")
    if RANDOM_SEED < 0:
        raise ValueError("RANDOM_SEED는 0 이상이어야 합니다.")
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
        str(TOTAL_COUNT),
        "--seed",
        str(RANDOM_SEED),
        "--outdir",
        str(output_dir),
    ]
    if not MAKE_PNG:
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
            str(TOTAL_COUNT),
        ],
        check=True,
    )

    print("\n완료")
    print(f"GDS directory: {output_dir}")


if __name__ == "__main__":
    main()
