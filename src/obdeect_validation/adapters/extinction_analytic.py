"""Beer–Lambert control; this case does not compare against sim_telarray."""

import math
from pathlib import Path

from obdeect_validation.adapters import common


def main():
    from obdeect.propagation import TransmissionTable

    # All sampled emission heights lie on this constant optical-depth plateau.
    path = Path("constant-optical-depth.dat")
    path.write_text("# H2= 2.000, H1= 3.000 4.000\n400 0.5 0.5\n")
    table = TransmissionTable(str(path), 2000.0)
    samples = []
    for zenith_deg in (0, 20, 40, 60):
        airmass = 1 / math.cos(math.radians(zenith_deg))
        expected = math.exp(-0.5 * airmass)
        for height_m in (3000, 3500, 4000, 10000):
            actual = float(table.transmission(400, height_m, airmass))
            if not math.isfinite(actual):
                raise ValueError("non-finite transmission")
            samples.append(
                {
                    "zenith_deg": zenith_deg,
                    "height_m": height_m,
                    "reference": expected,
                    "candidate": actual,
                }
            )
    common.write_json("samples.json", samples)
    _plot(samples, Path("extinction-analytic.png"))
    common.write_results(
        {
            "sample_count": len(samples),
            "max_transmission_absolute_error": max(
                abs(s["candidate"] - s["reference"]) for s in samples
            ),
        }
    )


def _plot(samples, path):
    """Transmission against emission height, one colour per zenith angle: analytic line, obdeect markers."""
    from obdeect_validation.adapters import plotting

    fig, axes = plotting.figure(1, 2, (10, 3.8))
    for i, zenith in enumerate(sorted({s["zenith_deg"] for s in samples})):
        rows = sorted(
            (s for s in samples if s["zenith_deg"] == zenith),
            key=lambda s: s["height_m"],
        )
        color = plotting.CYCLE[i % len(plotting.CYCLE)]
        heights = [s["height_m"] / 1000 for s in rows]
        axes[0, 0].plot(
            heights,
            [s["reference"] for s in rows],
            color=color,
            label=f"{zenith:g} deg",
        )
        axes[0, 0].plot(
            heights,
            [s["candidate"] for s in rows],
            "o",
            color=color,
            markerfacecolor="none",
        )
        axes[0, 1].plot(
            heights, [s["candidate"] - s["reference"] for s in rows], "o-", color=color
        )
    axes[0, 0].set(xlabel="emission height [km]", ylabel="transmission")
    axes[0, 0].legend(title="zenith (line: exp(-tau X), circles: obdeect)")
    axes[0, 1].set(xlabel="emission height [km]", ylabel="obdeect - analytic")
    plotting.save(fig, path)


if __name__ == "__main__":
    main()
