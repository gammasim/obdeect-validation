"""Command line interface for the validation runner."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .runner import ConfigError, brief, doctor, load_suite, report, run


def _json(value: object) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(prog="obdeect-validate")
    commands = result.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="list cases in a suite")
    listing.add_argument("--suite", required=True, type=Path)
    listing.add_argument("--json", action="store_true", dest="as_json")
    health = commands.add_parser("doctor", help="check the validation environment")
    health.add_argument("--environment", required=True, type=Path)
    health.add_argument("--json", action="store_true", dest="as_json")
    execution = commands.add_parser("run", help="run a suite or one case")
    choice = execution.add_mutually_exclusive_group(required=True)
    choice.add_argument("--suite", type=Path)
    choice.add_argument("--case", type=Path)
    execution.add_argument("--output", required=True, type=Path)
    execution.add_argument("--environment", type=Path)
    execution.add_argument("--non-interactive", action="store_true")
    execution.add_argument("--resume", action="store_true")
    for name in ("compare", "report"):
        command = commands.add_parser(name, help="show run results")
        command.add_argument("--run", required=True, type=Path)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "list":
            suite = load_suite(args.suite)
            cases = [{"id": c["id"], "name": c.get("name", c["id"])} for c in suite["cases"]]
            if args.as_json:
                _json({"suite": suite["id"], "cases": cases})
            else:
                for case in cases:
                    print(f"{case['id']}\t{case['name']}")
            return 0
        if args.command == "doctor":
            result = doctor(args.environment)
            if args.as_json:
                _json(result)
            else:
                for item in result["checks"]:
                    print(f"{item['status']}\t{item['name']}: {item['detail']}")
            return 0 if result["status"] == "PASS" else 3
        if args.command == "run":
            result = run(args.suite, args.case, args.output, args.environment, args.resume)
            _json(brief(result, args.output))
            if result.get("interrupted"):
                return 130
            return {"PASS": 0, "FAIL": 1, "ERROR": 2, "BLOCKED": 3}[result["status"]]
        result = report(args.run)
        _json(brief(result, args.run))
        return {"PASS": 0, "FAIL": 1, "ERROR": 2, "BLOCKED": 3}.get(result["status"], 2)
    except (ConfigError, OSError, ValueError) as exc:
        print(f"obdeect-validate: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
