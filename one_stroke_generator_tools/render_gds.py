"""GDS의 실제 금속 셀을 한 장의 PNG로 표시한다. 선 중심선을 다시 그리지 않는다."""

from pathlib import Path
import json

import gdstk
import numpy as np
from PIL import Image, ImageDraw, ImageFont, PngImagePlugin

PALETTE = np.array([[255, 255, 255], [209, 43, 47], [23, 104, 197], [121, 57, 158]], dtype=np.uint8)
LEFT, TOP, DISPLAY = 80, 100, 600


def read_masks(path, size_um):
    n = int(size_um / 5)
    masks = np.zeros((2, n, n), dtype=bool)
    lib = gdstk.read_gds(path)
    top = lib.top_level()
    if len(top) != 1 or top[0].references or top[0].paths:
        raise ValueError("flat rectangle GDS가 필요합니다.")
    for p in top[0].polygons:
        layer = {(39, 60): 0, (38, 40): 1}.get((p.layer, p.datatype))
        if layer is None:
            continue
        lo, hi = np.asarray(p.bounding_box())
        coords = np.concatenate((lo, hi)) / 5
        if len(p.points) != 4 or not np.allclose(coords, np.rint(coords), atol=1e-7, rtol=0):
            raise ValueError("금속이 5 um 직사각형 격자에 있지 않습니다.")
        x0, y0, x1, y1 = np.rint(coords).astype(int)
        if not (0 <= x0 < x1 <= n and 0 <= y0 < y1 <= n):
            raise ValueError("금속이 canvas를 벗어났습니다.")
        masks[layer, y0:y1, x0:x1] = True
    return masks, top[0].labels


def metal_rgb(masks):
    return PALETTE[(masks[0].astype(np.uint8) + 2 * masks[1].astype(np.uint8))[::-1]]


def render_gds(gds_path, png_path, size_um, subtitle=""):
    masks, labels = read_masks(gds_path, size_um)
    n = masks.shape[1]
    block = DISPLAY // n
    pixels = np.repeat(np.repeat(metal_rgb(masks), block, axis=0), block, axis=1)
    image = Image.new("RGB", (780, 785), "white")
    image.paste(Image.fromarray(pixels), (LEFT, TOP))
    draw = ImageDraw.Draw(image)
    # 기본 글꼴은 외부 폰트 설치 없이도 사용할 수 있다.
    import matplotlib
    font_path = Path(matplotlib.get_data_path()) / "fonts/ttf/DejaVuSans.ttf"
    font = ImageFont.truetype(str(font_path), 15)
    draw.text((LEFT, 15), f"{Path(gds_path).stem} | {size_um:g} x {size_um:g} um", fill="black", font=font)
    draw.text((LEFT, 37), subtitle, fill="black", font=font)
    draw.text((LEFT, 59), "M9 red | M8 blue | projected overlap purple", fill="black", font=font)
    draw.text((LEFT, 78), "5 um cells; actual metal/pads; via cuts omitted in PNG", fill="black", font=font)
    draw.rectangle((LEFT - 1, TOP - 1, LEFT + DISPLAY, TOP + DISPLAY), outline="black")
    for value in range(0, int(size_um) + 1, 50):
        pos = int(DISPLAY * value / size_um)
        draw.text((LEFT + pos - 5, TOP + DISPLAY + 8), str(value), fill="black", font=font)
        draw.text((LEFT - 30, TOP + DISPLAY - pos - 4), str(value), fill="black", font=font)
    for label in labels:
        # 포트 이름은 그림 밖에 적어 금속을 가리지 않는다.
        ycenter = label.origin[1] + 2.5
        y = TOP + DISPLAY - round(DISPLAY * ycenter / size_um)
        x = 2 if label.text.startswith("IN") else LEFT + DISPLAY + 10
        draw.text((x, y), label.text, fill="black", font=font)
    draw.text((LEFT, 755), "Layout preview only; no EM performance implied.", fill="black", font=font)
    info = PngImagePlugin.PngInfo()
    info.add_text("geometry_roi", json.dumps([LEFT, TOP, DISPLAY, DISPLAY]))
    info.add_text("size_um", str(size_um))
    info.add_text("grid_um", "5")
    image.save(png_path, pnginfo=info)


def preview_matches(gds_path, png_path, size_um):
    """PNG 금속 영역 전체를 GDS에서 복원한 RGB와 픽셀 단위로 대조한다."""
    masks, _ = read_masks(gds_path, size_um)
    block = DISPLAY // masks.shape[1]
    expected = np.repeat(np.repeat(metal_rgb(masks), block, axis=0), block, axis=1)
    with Image.open(png_path) as im:
        actual = np.asarray(im.convert("RGB"))[TOP:TOP + DISPLAY, LEFT:LEFT + DISPLAY]
    return np.array_equal(actual, expected)
