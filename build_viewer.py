"""Recalculate the animation data and build the interactive browser viewer.

Run this file from Spyder or with: python build_viewer.py
The force calculations come from the supplied part2a.py, not a second model.
"""

from pathlib import Path
import json
import re

import numpy as np
from plotly.offline import get_plotlyjs

from part2a import jacket_forces, PERIOD


ROOT = Path(__file__).resolve().parent
VIEWER = ROOT / "viewer"
DATA_FILE = ROOT / "generated" / "load_cycle.npz"


def sample_cycle():
    headings = np.arange(0.0, 360.0, 12.0)
    phases = np.linspace(0.0, 360.0, 100, endpoint=False)
    loads = np.zeros((len(phases), len(headings), len(headings)))

    for i, wave in enumerate(headings):
        print(f"Calculating wave direction {i + 1}/{len(headings)}...")
        for j, current in enumerate(headings):
            result = jacket_forces(wave, current, phases, 129)
            loads[:, i, j] = result["resultant"] / 1_000_000

    DATA_FILE.parent.mkdir(exist_ok=True)
    np.savez_compressed(
        DATA_FILE, headings=headings, phases=phases,
        load_MN=loads, period=PERIOD
    )
    return headings, phases, loads


def build_viewer():
    headings, phases, loads = sample_cycle()
    data = {
        "headings": headings.tolist(),
        "phases": phases.tolist(),
        "period": PERIOD,
        "loads": np.round(loads, 6).tolist(),
    }
    (VIEWER / "loads.js").write_text(
        "window.JACKET_LOAD_DATA="
        + json.dumps(data, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )
    (VIEWER / "plotly.min.js").write_text(get_plotlyjs(), encoding="utf-8")

    # Bundle the chart, scripts and styles into one portable HTML file.
    html = (VIEWER / "index.html").read_text(encoding="utf-8")
    css = (VIEWER / "styles.css").read_text(encoding="utf-8")
    html = html.replace(
        '<link rel="stylesheet" href="styles.css">',
        "<style>" + css + "</style>",
    )
    for name in ["plotly.min.js", "loads.js", "viewer.js"]:
        html = html.replace('<script defer src="' + name + '"></script>', "")
        script = (VIEWER / name).read_text(encoding="utf-8")
        script = script.replace("</script", "<\\/script")
        html = html.replace("</body>", "<script>" + script + "</script>\n</body>")
    html = re.sub(r'<nav aria-label="Downloads">.*?</nav>', "", html, flags=re.S)
    output = VIEWER / "jacket_load_interactive.html"
    output.write_text(html, encoding="utf-8")
    print(f"\nOpen this file in your browser: {output}")
    print("Drag to rotate, use the phase slider, or press Play cycle.")


if __name__ == "__main__":
    build_viewer()
