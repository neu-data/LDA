"""
book/render_diagrams.py
Draws each Mermaid diagram of the book as a PNG image.

    python book/render_diagrams.py en

Reads book/<lang>/mermaid/src/<id>.mmd (written by make_book_sources.R),
draws each diagram in its own page with the Mermaid library bundled with
Quarto, screenshots it with headless Microsoft Edge (or Chrome) at high
resolution and trims the surrounding white space. One page per diagram keeps
the mapping between source and image certain, and the large window avoids the
cropping seen with Quarto's own PNG export.

Requires: Python 3 with Pillow; Edge or Chrome; Quarto (for mermaid.min.js).
Set BROWSER_PATH or MERMAID_JS to override the default locations.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from PIL import Image, ImageChops

ROOT = Path(__file__).resolve().parent
LANG = sys.argv[1] if len(sys.argv) > 1 else "en"
MM_DIR = ROOT / LANG / "mermaid"
SRC_DIR = MM_DIR / "src"

BROWSER = os.environ.get(
    "BROWSER_PATH", r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")
MERMAID_JS = Path(os.environ.get(
    "MERMAID_JS",
    r"C:\Program Files\RStudio\resources\app\bin\quarto\share\formats\html\mermaid\mermaid.min.js"))

PAGE = """<!doctype html><html><head><meta charset="utf-8">
<style>body{{margin:0;padding:24px;background:#fff;font-family:"Segoe UI",Arial,sans-serif}}</style>
<script>{js}</script></head><body>
<pre class="mermaid">{code}</pre>
<script>mermaid.initialize({{startOnLoad:true, theme:"neutral",
  flowchart:{{useMaxWidth:false, htmlLabels:true}}}});</script>
</body></html>"""


def trim(path, pad=24):
    """Crop to the non-white content plus a small margin."""
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    box = ImageChops.difference(im, bg).getbbox()
    if box:
        l, t, r, b = box
        im = im.crop((max(l - pad, 0), max(t - pad, 0),
                      min(r + pad, im.width), min(b + pad, im.height)))
    im.save(path, optimize=True)


def render(mmd, png):
    code = mmd.read_text(encoding="utf-8").replace("&", "&amp;").replace("<", "&lt;")
    code = code.replace("&lt;br/>", "<br/>")          # keep line breaks in labels
    js = MERMAID_JS.read_text(encoding="utf-8")
    tmp = Path(tempfile.mkdtemp(prefix="mmd_"))
    page = tmp / "diagram.html"
    page.write_text(PAGE.format(js=js, code=code), encoding="utf-8")
    if png.exists():
        png.unlink()
    proc = subprocess.Popen(
        [BROWSER, "--headless", "--disable-gpu", f"--user-data-dir={tmp / 'profile'}",
         "--hide-scrollbars", "--force-device-scale-factor=2",
         "--window-size=1600,2200", "--virtual-time-budget=8000",
         f"--screenshot={png}", page.as_uri()],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    start = time.time()
    while time.time() - start < 90:                 # wait for the file to settle
        time.sleep(2)
        if png.exists() and png.stat().st_size > 0:
            size = png.stat().st_size
            time.sleep(2)
            if png.stat().st_size == size:
                break
    try:
        proc.wait(timeout=20)
    except subprocess.TimeoutExpired:
        proc.kill()
    shutil.rmtree(tmp, ignore_errors=True)
    if not png.exists():
        raise RuntimeError(f"No image produced for {mmd.name}")
    trim(png)


def main():
    sources = sorted(SRC_DIR.glob("*.mmd"))
    for mmd in sources:
        png = MM_DIR / (mmd.stem + ".png")
        render(mmd, png)
        print(f"{png.name}: {Image.open(png).size[0]} x {Image.open(png).size[1]} px")
    print(f"Rendered {len(sources)} diagrams to {MM_DIR}")


if __name__ == "__main__":
    main()
