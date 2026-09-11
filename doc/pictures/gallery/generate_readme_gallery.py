#!/usr/bin/env python3
"""
Regenerate the animated README gallery GIF from doc/pictures/gallery/.

Every image placed in doc/pictures/gallery/ is automatically picked up (sorted
by filename) and appended to doc/pictures/readme_gallery.gif, which is
embedded in the README's Gallery section. A static image becomes one held
frame; an already-animated source has all of its
frames spliced in at their original per-frame durations, so its motion is
preserved rather than collapsed to a single still. Run this script again
after adding, removing, or replacing images in that folder.

Usage: python3 doc/pictures/gallery/generate_readme_gallery.py
"""

from pathlib import Path

from PIL import Image, ImageSequence

GALLERY_DIR = Path(__file__).resolve().parent
REPO_ROOT = GALLERY_DIR.parent.parent.parent
SOURCE_DIR = GALLERY_DIR
OUTPUT_PATH = REPO_ROOT / "doc" / "pictures" / "readme_gallery.gif"

FRAME_WIDTH = 800
STATIC_HOLD_MS = 3000
# Many animated WebP encoders don't round-trip per-frame duration through
# Pillow; fall back to a plausible video rate instead of the (much longer)
# static-slide hold time.
ANIMATED_FALLBACK_DURATION_MS = 100
# Thin out already-animated sources to keep the merged GIF a reasonable size;
# duration is scaled up to match, so playback speed is unchanged.
ANIMATED_FRAME_STRIDE = 2
IMAGE_EXTENSIONS = {".png", ".jpg", ".jpeg", ".webp", ".gif", ".bmp"}


def frames_for(path):
    """Yield (frame, duration_ms) for every frame of an image, resized and
    palettized with its own local color table so unrelated slides don't
    smear into a shared global palette."""
    with Image.open(path) as im:
        aspect = im.height / im.width
        size = (FRAME_WIDTH, round(FRAME_WIDTH * aspect))
        is_animated = getattr(im, "n_frames", 1) > 1
        stride = ANIMATED_FRAME_STRIDE if is_animated else 1
        fallback_duration = ANIMATED_FALLBACK_DURATION_MS if is_animated else STATIC_HOLD_MS

        for index, frame in enumerate(ImageSequence.Iterator(im)):
            if index % stride != 0:
                continue
            duration = frame.info.get("duration") or fallback_duration
            # convert("RGB") alone drops the alpha channel without flattening
            # it, leaving whatever (undefined) RGB value sat behind
            # transparent pixels - which differs frame to frame and shows up
            # as flicker. Composite onto opaque white first instead.
            rgba = frame.convert("RGBA")
            on_white = Image.alpha_composite(Image.new("RGBA", rgba.size, "white"), rgba)
            resized = on_white.convert("RGB").resize(size, Image.LANCZOS)
            palettized = resized.convert("P", palette=Image.ADAPTIVE, colors=256)
            yield palettized, duration * stride


def load_frames():
    paths = sorted(p for p in SOURCE_DIR.iterdir() if p.suffix.lower() in IMAGE_EXTENSIONS)
    if not paths:
        raise SystemExit(f"No images found in {SOURCE_DIR}")

    frames, durations = [], []
    for path in paths:
        for frame, duration in frames_for(path):
            frames.append(frame)
            durations.append(duration)
    return frames, durations


def main():
    frames, durations = load_frames()
    frames[0].save(
        OUTPUT_PATH,
        save_all=True,
        append_images=frames[1:],
        duration=durations,
        loop=0,
        optimize=True,
        # Every frame is a full, opaque replacement of the canvas (no
        # sprite-style partial updates), so disposal must leave the previous
        # frame in place rather than clearing to the GIF background color -
        # clearing between frames flashes a wrong color from each frame's
        # local palette and looks like a flickering/transparent background.
        disposal=1,
    )
    print(f"Wrote {OUTPUT_PATH} from {len(frames)} frame(s)")


if __name__ == "__main__":
    main()
