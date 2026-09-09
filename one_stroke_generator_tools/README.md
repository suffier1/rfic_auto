# One-stroke GDS Generator

M9·M8 한붓그리기 GDS를 생성하고 구조·연결 상태를 검사한다.
EM 시뮬레이션은 실행하지 않는다. 생성한 GDS를 기존 연구실 EM 작업에 입력한다.

## 파일 구성

- `easy_generate_one_stroke.py`: 생성 설정을 읽어 GDS 생성과 검증을 순서대로 실행한다.
- `generate_one_stroke_gds.py`: one-stroke 형상과 포트 배치를 생성한다.
- `em_gds_contract.py`: GDS 레이어, 포트, VIA 규격을 정의하고 GDS 파일을 기록한다.
- `verify_one_stroke_gds.py`: GDS 구조, 포트, VIA, 연결 및 절연 상태를 검사한다.
- `analyze_sparams.py`: manifest와 s4p를 숫자 ID로 연결해 지정 주파수의 차동 성능을 집계한다.
- `test_one_stroke.py`: 기존 형상 재현·분할 생성·분석 회귀 검사다.
- `requirements_one_stroke.txt`: 필요한 Python 패키지 목록이다.
- `EM_WORKFLOW.md`: 77 GHz 분석 근거, 생성 방향, EM 실행·반환 조건이다.

## 설치

Python 3.10 이상을 사용한다. 이 디렉터리에서 실행한다.

```bash
python -m pip install -r requirements_one_stroke.txt
```

## 먼저 30개 생성·검사

```bash
python easy_generate_one_stroke.py --n 30 --seed 42 --family serpentine --outdir ../batches/serpentine_check30
```

생성 후 전수검사가 자동 실행된다. PNG는 기본적으로 생략하며, 필요하면 `--png`를 추가한다.
기존 출력 디렉터리는 덮어쓰지 않는다.

## 1만 개 생성과 분할 실행

```bash
python easy_generate_one_stroke.py --n 10000 --seed 42 --family serpentine --start-index 0 --outdir ../batches/serpentine_seed42_00000_09999
```

다음 1만 개는 같은 seed를 유지하고 시작 ID를 늘린다.

```bash
python easy_generate_one_stroke.py --n 10000 --seed 42 --family serpentine --start-index 10000 --outdir ../batches/serpentine_seed42_10000_19999
```

시작 ID를 `0, 10000, …, 90000`으로 바꾸면 총 10만 개를 10개 batch로 나눌 수 있다.
같은 seed·family·index는 같은 형상을 생성하므로 ID 구간을 겹치지 않게 한다.
분할 여부는 형상에 영향을 주지 않는다. 다른 seed 또는 family의 같은 ID는 다른 데이터이므로 디렉터리를 섞지 않는다.

기존 세 형상을 고르게 생성하려면 `--family all`을 사용한다.
`--family large_rect`, `--family deep_loop`도 가능하다.

## 설정값을 직접 수정하는 경우

`easy_generate_one_stroke.py` 상단의 값을 수정한다.

```python
TOTAL_COUNT = 1000
RANDOM_SEED = 8
OUTPUT_FOLDER_NAME = "one_stroke_gds_seed_8_new"
MAKE_PNG = False
FAMILY = "all"
START_INDEX = 0
```

- `TOTAL_COUNT`: 생성 개수
- `RANDOM_SEED`: 난수 seed
- `OUTPUT_FOLDER_NAME`: 출력 디렉터리 이름
- `MAKE_PNG`: PNG 미리보기 생성 여부
- `FAMILY`: `all / large_rect / deep_loop / serpentine`
- `START_INDEX`: 전역 ID 시작값. 같은 seed로 분할 생성할 때 이전 구간 다음 번호

명령행 인자를 함께 주면 명령행 값이 우선한다.
`OUTPUT_FOLDER_NAME`은 코드 디렉터리 기준이고, 명령행의 `--outdir`은 현재 작업 디렉터리 기준이다.

## 실행

```bash
python easy_generate_one_stroke.py
```

## 출력 파일

| 파일 | 내용 |
|---|---|
| `difftx_XXXXX.gds` | EM 입력용 GDS |
| `manifest.csv` | seed, 전역 ID, 형상·배치, 길이·좌표, serpentine PRI/SEC 횡단 구간 수 |
| `dataset_meta.json` | 생성 버전·설정·소스 해시와 완료 상태 |
| `verification.json` | GDS 규격·연결 전수검사 결과 |
| `difftx_XXXXX.png` | `--png`를 지정한 경우의 미리보기 |

`dataset_meta.json`의 `status=complete`와 `verification.json`의 `passed=true`를 확인한다.
이는 생성과 구조 검사의 완료이며 EM 검증 완료를 뜻하지 않는다.
중단된 batch는 `status=generating`으로 남는다. 자동 resume/압축/삭제는 하지 않으며,
중단된 batch를 EM에 넣지 말고 새 출력 디렉터리에서 같은 설정으로 다시 생성한다.

## EM 결과 분석

EM 조건과 반환 파일은 [EM_WORKFLOW.md](EM_WORKFLOW.md)를 따른다.

```bash
python analyze_sparams.py --manifest ../batches/serpentine_seed42_00000_09999/manifest.csv --sparams ../em_results/serpentine_seed42_00000_09999 --freq-ghz 77 --outdir ../analysis_results/serpentine_seed42_00000_09999
```

`--sparams`에는 s4p 디렉터리 또는 tar.gz를 지정한다. 압축을 풀지 않고 읽는다.
Touchstone 1의 4-port, 포트별 50 Ω 기준을 지원한다. 포트 1·2가 PRI, 3·4가 SEC여야 한다.
출력은 샘플별 `per_sample.csv`와 형상·배치별 `summary.json`이다.

- `target_source=direct_em`: 해당 주파수를 직접 EM으로 계산한 값.
- `target_source=complex_linear_interpolation`: 해당 지점이 없어 양옆의 복소 S를 보간한 값. 직접 EM 값이 아님.
- `*_pass_direct`: 직접 EM 지점의 통과 여부. 보간이면 빈칸.
- `*_pass_interpolated`: 보간값의 통과 여부. 직접 EM이면 빈칸.
- `*_pass_both_brackets`: 양옆 샘플 지점이 모두 통과하는지. 연속 대역 전체 통과를 뜻하지 않음.
- `missing_indices`: manifest에 있으나 s4p가 없는 ID. 통과율 분모는 반환된 s4p 수로 확인.

좋은 결과만 남겨 원본 s4p를 버리지 않는다. 조건부 모델 학습에는 정상적인 낮은 성능 사례도 사용한다.

## 형상 변경과 검사

형상을 변경하려면 `generate_one_stroke_gds.py`의 `build_*_route()` 함수를
수정한다. 레이어, 포트 또는 VIA 규격을 변경할 때만 `em_gds_contract.py`를
수정한다.

기본 경로 분포는 v3 그대로다. `serpentine`은 권선마다 4/6 횡단 구간을 독립적으로 선택하고,
세 배치 방식을 고르게 생성한다. 성능에 따른 거절 선별은 생성기에 넣지 않았다.
`--family all --start-index 0`은 기존 seed/index의 형상을 재현한다.
기본 GDS 검사는 foundry DRC signoff나 EM 정확성 검사가 아니다.

```bash
python -m unittest -v test_one_stroke
```
