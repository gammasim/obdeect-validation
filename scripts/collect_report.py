"""Copy the small, reviewable part of a validation run into reports/.

Keeps verdicts, metrics, details and figures; drops photon files, simulator
output, logs and generated model configurations.

    python scripts/collect_report.py results/full-integrated reports/2026-11-13
"""

import json
import shutil
import sys
from pathlib import Path

KEEP_NAMES = {"metrics.json", "details.json"}
KEEP_SUFFIXES = {".png", ".svg"}
TOP_LEVEL = ("summary.json", "junit.xml", "report.html")
MAX_BYTES = 2_000_000


def collect(run: Path, report: Path) -> list[Path]:
    copied = []
    report.mkdir(parents=True, exist_ok=True)
    for name in TOP_LEVEL:
        if (run / name).is_file():
            shutil.copy2(run / name, report / name)
            copied.append(report / name)
    for path in sorted(run.rglob("*")):
        relative = path.relative_to(run)
        if (
            not path.is_file()
            or "cfg" in relative.parts
            or path.stat().st_size > MAX_BYTES
        ):
            continue
        if path.name in KEEP_NAMES or path.suffix in KEEP_SUFFIXES:
            target = report / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            copied.append(target)
    return copied


def main() -> None:
    run, report = Path(sys.argv[1]), Path(sys.argv[2])
    copied = collect(run, report)
    total = sum(path.stat().st_size for path in copied)
    print(json.dumps({"files": len(copied), "bytes": total}))


if __name__ == "__main__":
    main()
