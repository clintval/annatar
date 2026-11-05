"""Render the README figure by running the README's own MuData code under the figure style.

Run from the repository root, after the README's `annatar` command, in an environment with muon.
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

REPO = Path(__file__).resolve().parents[2]
STYLE = {
    "font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8,
    "text.color": "#0b0b0b",
    "axes.titlecolor": "#0b0b0b",
    "axes.titlelocation": "left",
    "axes.titlepad": 8,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "figure.dpi": 200,
    "figure.figsize": (3.0, 2.9),
    "legend.fontsize": 7,
}


def main() -> None:
    """Write the figure the README's code draws."""
    readme = (REPO / "README.md").read_text()
    section = readme.split("## MuData Example\n")[1].split("\n## ")[0]
    code = re.search(r"```python\n(.*?)```", section, flags=re.DOTALL)
    plt.switch_backend("Agg")
    plt.rcParams.update(STYLE)
    sc.settings.autoshow = False
    exec(code.group(1), {})
    plt.gcf().savefig(
        REPO / ".github" / "img" / "annatar_example.png",
        dpi=200,
        bbox_inches="tight",
        pad_inches=0.05,
        facecolor="white",
    )


if __name__ == "__main__":
    main()
