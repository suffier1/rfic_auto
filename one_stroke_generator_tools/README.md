# One-stroke GDS Generator

One-stroke 형상의 GDS를 생성하고, 생성된 파일의 구조와 연결 상태를 검사한다.

## 파일 구성

- `easy_generate_one_stroke.py`: 생성 설정을 읽어 GDS 생성과 검증을 순서대로 실행한다.
- `generate_one_stroke_gds.py`: one-stroke 형상과 포트 배치를 생성한다.
- `em_gds_contract.py`: GDS 레이어, 포트, VIA 규격을 정의하고 GDS 파일을 기록한다.
- `verify_one_stroke_gds.py`: GDS 구조, 포트, VIA, 연결 및 절연 상태를 검사한다.
- `requirements_one_stroke.txt`: 필요한 Python 패키지 목록이다.

## 설치

```bash
python -m pip install -r requirements_one_stroke.txt
```

## 설정

`easy_generate_one_stroke.py` 상단의 값을 수정한다.

```python
TOTAL_COUNT = 1000
RANDOM_SEED = 8
OUTPUT_FOLDER_NAME = "one_stroke_gds_seed_8_new"
MAKE_PNG = False
```

- `TOTAL_COUNT`: 생성 개수
- `RANDOM_SEED`: 난수 seed
- `OUTPUT_FOLDER_NAME`: 출력 디렉터리 이름
- `MAKE_PNG`: PNG 미리보기 생성 여부

## 실행

```bash
python easy_generate_one_stroke.py
```

지정한 디렉터리에 GDS와 `manifest.csv`가 생성된다. 생성이 끝나면 모든 GDS를
자동으로 검증한다. 기존 출력 디렉터리는 덮어쓰지 않는다.

형상을 변경하려면 `generate_one_stroke_gds.py`의 `build_*_route()` 함수를
수정한다. 레이어, 포트 또는 VIA 규격을 변경할 때만 `em_gds_contract.py`를
수정한다.
