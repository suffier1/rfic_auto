# Batch GDS Generator

두 권선의 크기·경로를 제어해 one-stroke GDS 배치를 만든다. 경로 부품과 GDS 규격은 `../one_stroke_generator_tools`를 그대로 쓴다.

- `gen_asym.py`: 1차/2차 목표 면적비를 [1/3.5, 3.5]에서 로그 균등으로 뽑아 두 권선 크기를 정한다.
- `gen_mix.py`: 크기가 비슷한 3,000개(형상 3종 × 크기 3단계 × 상대 배치 3종)와 비대칭 2,000개. 굴곡 깊이·행 수·행 간격·반환 위치를 따로 뽑는다.
- `generation_gates.py`: 생성 중에 거르는 검사(한 에지에 pin 하나, 벽 단자, 구멍). 떨어진 경로는 수리하지 않고 다시 뽑는다.

## 실행

이 디렉터리에서 실행한다. 기존 출력 폴더는 덮어쓰지 않는다.

```bash
python -m pip install -r requirements.txt
python gen_asym.py --n 5000 --seed 2026 --outdir ../batches/asym5000
python gen_mix.py --seed 2027 --outdir ../batches/mix5000
python gen_mix.py --pilot --seed 2027 --outdir ../batches/mix_pilot    # 구성별 대표 78개
```

출력: `<outdir>/gds/*.gds`, `manifest.csv`(샘플별 설정), `gen_meta.json`.
같은 seed와 코드면 같은 형상이 나온다. GDS 파일에는 저장 시각이 들어가므로 파일 해시는 매번 다르다.
