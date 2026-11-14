"""Figures shared by the adapters: one palette, one style, helpers for the common panels.

Colours carry meaning and are the same in every plot (Okabe-Ito, colour-blind safe):
obdeect blue, the reference program vermilion, analytic expectations black. Cases and
telescopes get their own fixed colours so a case is recognisable across figures.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

NATIVE = "#0072B2"
REFERENCE = "#D55E00"
ANALYTIC = "#000000"
LIMIT = "#999999"

# Okabe-Ito colours without black and yellow (poor contrast on white).
CYCLE = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#6A3D9A")

# Fixed colour per case family and per telescope type, in suite order.
CASE_COLORS = {
    "stream": "#0072B2",
    "propagation": "#009E73",
    "sources": "#E69F00",
    "optics": "#CC79A7",
    "chain": "#D55E00",
    "performance": "#56B4E9",
    "lst": "#0072B2",
    "mst": "#009E73",
    "sst": "#D55E00",
}

plt.rcParams.update(
    {
        "axes.prop_cycle": matplotlib.cycler(color=CYCLE),
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "legend.frameon": False,
        "legend.fontsize": "small",
        "figure.dpi": 100,
        "savefig.dpi": 130,
        "savefig.bbox": "tight",
    }
)


def case_color(name: str) -> str:
    """Colour of a case family or telescope ('optics', 'LST', ...); stable for unknown names."""
    key = name.lower()
    for known, color in CASE_COLORS.items():
        if known in key:
            return color
    return CYCLE[sum(map(ord, key)) % len(CYCLE)]


def figure(rows: int = 1, columns: int = 1, size: tuple[float, float] | None = None):
    """Figure and axes (always a 2-d array of axes for rows x columns > 1)."""
    size = size or (4.6 * columns + 0.6, 3.4 * rows + 0.4)
    fig, axes = plt.subplots(rows, columns, figsize=size, squeeze=False)
    return fig, axes if rows * columns > 1 else axes[0, 0]


def save(fig, path) -> None:
    fig.tight_layout()
    fig.savefig(Path(path))
    plt.close(fig)


def message(path, text: str) -> None:
    """A figure that only states why there is nothing to plot."""
    fig, axis = figure(size=(7, 2))
    axis.axis("off")
    axis.text(0.02, 0.5, text, va="center", wrap=True)
    save(fig, path)


def histograms(
    axis,
    samples,
    xlabel: str,
    bins: int = 60,
    density: bool = True,
    clip: tuple[float, float] | None = (0.001, 0.999),
) -> None:
    """Overlaid step histograms on a common range.

    samples: (label, values, colour) tuples; empty samples are skipped.
    """
    samples = [(label, np.asarray(v), c) for label, v, c in samples if len(v)]
    if samples:
        everything = np.concatenate([v for _, v, _ in samples])
        low, high = (
            np.quantile(everything, clip)
            if clip
            else (everything.min(), everything.max())
        )
        if not high > low:
            low, high = low - 0.5e-6, low + 0.5e-6
        for label, values, color in samples:
            axis.hist(
                values,
                bins=bins,
                range=(low, high),
                density=density,
                histtype="stepfilled",
                color=color,
                alpha=0.10,
                linewidth=0,
            )
            axis.hist(
                values,
                bins=bins,
                range=(low, high),
                density=density,
                histtype="step",
                color=color,
                linewidth=1.4,
                label=label,
            )
        axis.legend()
    axis.set(xlabel=xlabel, ylabel="density" if density else "entries")


def camera_bars(axis, native, reference, reference_label: str) -> None:
    """Mean p.e. per triggered event of the two samples with their standard errors."""
    labels, colors = ["obdeect", reference_label], [NATIVE, REFERENCE]
    means, errors = [], []
    for values in (native, reference):
        values = np.asarray(values, dtype=float)
        means.append(float(values.mean()) if values.size else 0.0)
        errors.append(
            float(values.std(ddof=1) / math.sqrt(values.size))
            if values.size > 1
            else 0.0
        )
    axis.bar(labels, means, yerr=errors, color=colors, capsize=4, alpha=0.85)
    axis.set(ylabel="mean p.e. per triggered event")
    axis.grid(axis="x", visible=False)
    for i, values in enumerate((native, reference)):
        if len(values) == 0:
            axis.text(
                i, 0, "no triggered events", ha="center", va="bottom", rotation=90
            )
