# Requires scanpy, muon and matplotlib, which the pixi environment does not provide.
"""Render the README figure by running the README's own MuData code under the figure style.

Run from the repository root, after the README's `annatar` command.
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

REPO = Path(__file__).resolve().parents[2]
DPI = 400
SPINE_COLOR = "#e1e0d9"
STYLE = {
    "font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8,
    "text.color": "#0b0b0b",
    "axes.titlecolor": "#0b0b0b",
    "axes.titlelocation": "left",
    "axes.titlepad": 8,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "figure.dpi": DPI,
    "figure.figsize": (3.0, 2.9),
    "legend.fontsize": 7,
}


def frame(axes: plt.Axes) -> None:
    """Draw gray left and bottom axes with no labels or ticks."""
    axes.set_axis_on()
    axes.set_frame_on(True)
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_xlabel("")
    axes.set_ylabel("")
    axes.spines[["top", "right"]].set_visible(False)
    axes.spines[["left", "bottom"]].set_visible(True)
    axes.spines[["left", "bottom"]].set_color(SPINE_COLOR)
    axes.spines[["left", "bottom"]].set_linewidth(0.8)


def main() -> None:
    """Write the figure the README's code draws."""
    readme = (REPO / "README.md").read_text()
    section = readme.split("## MuData Example\n")[1].split("\n## ")[0]
    code = re.search(r"```python\n(.*?)```", section, flags=re.DOTALL)
    plt.switch_backend("Agg")
    plt.rcParams.update(STYLE)
    sc.settings.autoshow = False
    exec(code.group(1), {})
    for axes in plt.gcf().axes:
        frame(axes)
    plt.gcf().savefig(
        REPO / ".github" / "img" / "annatar_example.png",
        dpi=DPI,
        bbox_inches="tight",
        pad_inches=0.15,
        facecolor="white",
    )


if __name__ == "__main__":
    main()
