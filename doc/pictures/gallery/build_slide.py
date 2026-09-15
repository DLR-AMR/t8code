#!/usr/bin/env python
"""Compose one new slide for the README gallery.

Layout: white canvas, artwork on the left, title set in the whitespace on
the right. This builds exactly one new item - it does not know about or
rebuild any existing slide.

Usage:
    python3 build_gallery.py <input-file> <output> "<title>" \
        [--subtitle TEXT] [--fps N] [--quality Q]

<input-file> is a still image for a static slide, or a video/GIF for an
animated one - pass --fps to treat it as animated and set the frame
extraction rate. <output> names the resulting <output>.webp.
Use "\\n" in <title> to force a line break, otherwise it wraps to fit.
--subtitle adds a smaller caption line below the title (e.g. for a
citation).

After running this, wire <output> into wherever the README references the
gallery (e.g. re-run generate_readme_gallery.py to refresh the merged GIF).
"""
import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

GALLERY_DIR = Path(__file__).resolve().parent

FFMPEG = os.environ.get("T8_GALLERY_FFMPEG", shutil.which("ffmpeg") or "ffmpeg")
# Liberation Sans is the standard metric-compatible free substitute where
# Arial itself isn't installed/licensed (e.g. package "fonts-liberation" on
# Debian/Ubuntu); override via the env var if your system keeps it elsewhere.
FONT = os.environ.get(
    "T8_GALLERY_FONT",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf")
SUBTITLE_FONT = os.environ.get(
    "T8_GALLERY_SUBTITLE_FONT",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf")
SUBTITLE_FG = (140, 140, 140)

# Rendered at 2x (1920x640) so slides stay sharp on HiDPI screens.
CANVAS_W, CANVAS_H = 1920, 640
IMG_BOX = (40, 32, 1200, 608)       # left artwork area
TITLE_BOX = (1232, 48, 1880, 592)   # right whitespace reserved for the title
BG = (255, 255, 255)
FG = (38, 38, 38)


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


def flatten(img):
    """Composite any alpha onto white."""
    img = img.convert("RGBA")
    bg = Image.new("RGBA", img.size, BG + (255,))
    return Image.alpha_composite(bg, img).convert("RGB")


def wrap(text, font, max_w, draw):
    """Honour explicit \n, then greedily wrap each paragraph."""
    lines = []
    for para in text.split("\n"):
        words, cur = para.split(), ""
        for w in words:
            trial = (cur + " " + w).strip()
            if draw.textlength(trial, font=font) <= max_w or not cur:
                cur = trial
            else:
                lines.append(cur)
                cur = w
        lines.append(cur)
    return lines


def draw_title(canvas, text, subtitle=None):
    draw = ImageDraw.Draw(canvas)
    x0, y0, x1, y1 = TITLE_BOX
    max_w, max_h = x1 - x0, y1 - y0

    reserved, gap, sub_font = 0, 0, None
    if subtitle:
        sub_font = ImageFont.truetype(SUBTITLE_FONT, 30)
        gap = 16
        reserved = gap + int(sub_font.size * 1.28)

    for size in range(80, 26, -1):
        font = ImageFont.truetype(FONT, size)
        lines = wrap(text, font, max_w, draw)
        lh = int(size * 1.28)
        if (all(draw.textlength(l, font=font) <= max_w for l in lines)
                and lh * len(lines) <= max_h - reserved):
            break
    total_h = lh * len(lines)
    y = y0 + (max_h - total_h - reserved) // 2
    cx = (x0 + x1) / 2
    for line in lines:
        w = draw.textlength(line, font=font)
        draw.text((cx - w / 2, y), line, font=font, fill=FG)
        y += lh
    if subtitle:
        y += gap
        w = draw.textlength(subtitle, font=sub_font)
        draw.text((cx - w / 2, y), subtitle, font=sub_font, fill=SUBTITLE_FG)


def compose(img, title, subtitle=None):
    """Fit the (already flattened) image into the left box, title on the right."""
    bx0, by0, bx1, by1 = IMG_BOX
    bw, bh = bx1 - bx0, by1 - by0
    scale = min(bw / img.width, bh / img.height)
    img = img.resize((max(1, round(img.width * scale)),
                      max(1, round(img.height * scale))), Image.LANCZOS)
    canvas = Image.new("RGB", (CANVAS_W, CANVAS_H), BG)
    canvas.paste(img, (bx0 + (bw - img.width) // 2,
                       by0 + (bh - img.height) // 2))
    draw_title(canvas, title, subtitle=subtitle)
    return canvas


def build_still(input_path, output, title, subtitle):
    img = flatten(Image.open(input_path))
    out = GALLERY_DIR / f"{output}.webp"
    # Lossless WebP: bit-exact like PNG but noticeably smaller.
    compose(img, title, subtitle).save(out, lossless=True, quality=100, method=6)
    return out


def build_anim(input_path, output, title, subtitle, fps, quality):
    work = Path(tempfile.mkdtemp(prefix="t8_gallery_"))
    raw_dir, cmp_dir = work / "raw", work / "cmp"
    raw_dir.mkdir()
    cmp_dir.mkdir()

    run([FFMPEG, "-y", "-loglevel", "error", "-i", str(input_path),
         "-vf", f"fps={fps}", str(raw_dir / "%04d.png")])
    names = sorted(os.listdir(raw_dir))

    for name in names:
        img = flatten(Image.open(raw_dir / name))
        compose(img, title, subtitle).save(cmp_dir / name)

    # Animated WebP - truecolour and far smaller than an equivalent GIF.
    # Encoded via Pillow so kmin/kmax are reachable: libwebp's default
    # inter-frame coding re-quantises only the changed rectangle, which
    # leaves visible seams wherever a smooth gradient is partly redrawn.
    # Forcing every frame to be a keyframe removes them.
    out = GALLERY_DIR / f"{output}.webp"
    imgs = [Image.open(cmp_dir / n).convert("RGB") for n in names]
    imgs[0].save(out, save_all=True, append_images=imgs[1:],
                 duration=round(1000 / fps), loop=0, quality=quality,
                 method=6, kmin=1, kmax=1, minimize_size=False)

    shutil.rmtree(work, ignore_errors=True)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="Source image, GIF, or video file.")
    parser.add_argument("output", help="Output slide name (produces <output>.webp).")
    parser.add_argument("title", help='Title text; use "\\n" to force a line break.')
    parser.add_argument("--subtitle", default=None,
                         help="Optional smaller caption drawn below the title.")
    parser.add_argument("--fps", type=int, default=None,
                         help="Treat input as animated and extract frames at this rate.")
    parser.add_argument("--quality", type=int, default=80,
                         help="WebP quality for an animated slide (stills are always lossless).")
    args = parser.parse_args()

    if args.fps:
        out = build_anim(args.input, args.output, args.title, args.subtitle,
                          args.fps, args.quality)
    else:
        out = build_still(args.input, args.output, args.title, args.subtitle)
    print(f"{out.name}: {out.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
