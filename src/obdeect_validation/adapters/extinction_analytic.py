"""Beer–Lambert control; this case does not compare against sim_telarray."""

import json
import math
from pathlib import Path


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
                dict(
                    zenith_deg=zenith_deg,
                    height_m=height_m,
                    reference=expected,
                    candidate=actual,
                )
            )
    Path("samples.json").write_text(json.dumps(samples, indent=2) + "\n")
    Path("metrics.json").write_text(
        json.dumps(
            {
                "sample_count": len(samples),
                "max_transmission_absolute_error": max(
                    abs(s["candidate"] - s["reference"]) for s in samples
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
