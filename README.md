# RFIC one-stroke 데이터 생성

두 금속층의 한붓그리기 GDS 생성, 구조 검사, EM 결과 집계 코드다.
EM 시뮬레이션 자체는 각 연구실의 기존 환경에서 실행한다.

```bash
git clone https://github.com/suffier1/rfic_auto.git
cd rfic_auto/one_stroke_generator_tools
python -m pip install -r requirements_one_stroke.txt
python easy_generate_one_stroke.py --n 30 --seed 42 --family serpentine --outdir ../batches/serpentine_check30
```

이미 받은 저장소를 갱신할 때는 저장소 안에서 `git pull --ff-only`를 실행한다.
직접 수정한 파일과 충돌하면 수정본을 보존하고 차이를 확인한다.

- [생성·분할 실행·설정·파일 역할](one_stroke_generator_tools/README.md)
- [77 GHz 분석 근거와 EM 실행·반환 안내](one_stroke_generator_tools/EM_WORKFLOW.md)

생성 GDS와 EM 결과는 Git에 포함하지 않는다. 압축·업로드도 자동 실행하지 않는다.
