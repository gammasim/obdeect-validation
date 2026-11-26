"""obdeect command lines and the standard air-shower chain setup.

Every adapter builds its obdeect command lines here instead of spelling out argument lists.

``obdeect(*command, **options)`` gives ``["obdeect", *command, --option value ...]``. Option names map
to flags: ``i`` and ``o`` become ``-i`` and ``-o``, any other name ``--name`` with ``_`` replaced by
``-`` (``site_altitude`` -> ``--site-altitude``). Values: ``True`` gives the bare flag, ``None`` and
``False`` omit the option, a list or tuple gives several values after one flag
(``wavelength_range=(300, 600)`` -> ``--wavelength-range 300 600``), anything else one value.
Options appear in keyword order, so a command line keeps the order it was written in. Values are not
converted to strings: ``common.command`` and ``endtoend_common.pipeline`` do that.

For an option name that is not a Python identifier (or must be repeated, such as ``--set``), pass
``extra=[...]``: those arguments are appended as they are.

``input_argument(parser)`` adds the ``--input`` option of the chain cases (CORSIKA IACT file,
default the 20-100 GeV gamma showers).

``prepare_chain`` makes the inputs every air-shower chain case starts from: the Cherenkov band of
the CORSIKA run, the transmission table cut below sim_telarray's MAX_LAMBDA, a one-telescope camera
map and the compiled LSTN-01 model.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from obdeect_validation.adapters import atmosphere_chain as atm
from obdeect_validation.adapters import common


def _flag(name: str) -> str:
    if name in ("i", "o"):
        return f"-{name}"
    return "--" + name.replace("_", "-")


def obdeect(*command: Any, extra=(), executable="obdeect", **options: Any) -> list:
    """[executable, *command, options..., *extra] (see the module docstring).

    executable: "obdeect" (found on PATH when run) or a resolved path such as common.tool("obdeect").
    """
    argv: list = [executable, *command]
    for name, value in options.items():
        if value is None or value is False:
            continue
        argv.append(_flag(name))
        if value is True:
            continue
        if isinstance(value, (list, tuple)):
            argv.extend(value)
        else:
            argv.append(value)
    argv.extend(extra)
    return argv


def compile_model(
    output,
    model: str = "LSTN-01",
    version: str = "7.0.0",
    *,
    seed=1,
    zenith_deg=None,
    simulation_models=None,
) -> list:
    """obdeect ray-tracing compile of model/version into output (zenith_deg omitted when None)."""
    return obdeect(
        "ray-tracing",
        "compile",
        simulation_models=(
            common.models_root() if simulation_models is None else simulation_models
        ),
        model=model,
        version=version,
        seed=seed,
        zenith_deg=zenith_deg,
        output=output,
    )


def from_corsika(source, **options) -> list:
    """obdeect stream from-corsika source (the stream goes to standard output)."""
    return obdeect("stream", "from-corsika", source, **options)


def propagation_run(i, o, **options) -> list:
    """obdeect propagation run with options, then -i i -o o."""
    return obdeect("propagation", "run", **options, i=i, o=o)


def weighted_propagation(cut, band, i, o, *, seed, **options) -> list:
    """The chain propagation: cut table, no in-telescope path, site altitude, weighted, band."""
    return propagation_run(
        i,
        o,
        transmission=cut,
        pathlen="0",
        site_altitude=atm.SITE_ALTITUDE_M,
        seed=seed,
        mode="weighted",
        wavelength_range=tuple(band),
        **options,
    )


def trace_stream(model, i, o, *, seed, **options) -> list:
    """obdeect ray-tracing trace-stream --model model --seed seed options..., then -i i -o o."""
    return obdeect(
        "ray-tracing", "trace-stream", model=model, seed=seed, **options, i=i, o=o
    )


def to_simtel(i, o, **options) -> list:
    """obdeect stream to-simtel i -o o options... (i "-" reads standard input)."""
    return obdeect("stream", "to-simtel", i, o=o, **options)


#: CORSIKA input of the air-shower chain cases.
SHOWERS_20_100 = common.VALIDATION / "inputs/gamma-20deg-20-100GeV.iact"


def input_argument(parser, default=SHOWERS_20_100) -> None:
    """Add --input (CORSIKA IACT file, a Path) to an argparse parser."""
    parser.add_argument("--input", type=Path, default=default)


@dataclass
class Chain:
    """Inputs of the standard air-shower chain (see prepare_chain)."""

    root: Path
    source: Path
    band: tuple[float, float]
    cut: Path
    camera_map: Path | None
    model: Path | None


def prepare_chain(
    root: Path,
    source: Path,
    *,
    model: str | None = "LSTN-01",
    version: str = "7.0.0",
    zenith_deg=None,
    compile_seed=1,
    camera_map_text: str | None = "1 1 0.0\n",
) -> Chain:
    """Band, cut table, camera map and (model not None) compiled model of the standard chain in root.

    The band comes from <source>.corsika.log; the cut table is common.TRANSMISSION_TABLE cut at
    atm.MAX_LAMBDA_NM (cut-transmission.dat); the camera map is camera-map.txt; the model is
    model.json, compiled with compile.log. camera_map_text None: no camera map (Chain.camera_map None).
    """
    root, source = Path(root), Path(source)
    band = atm.cherenkov_band(
        source.with_suffix(".corsika.log").read_text(errors="replace")
    )
    cut = root / "cut-transmission.dat"
    cut.write_text(
        atm.cut_table(common.TRANSMISSION_TABLE.read_text(), atm.MAX_LAMBDA_NM)
    )
    camera_map = None
    if camera_map_text is not None:
        camera_map = root / "camera-map.txt"
        camera_map.write_text(camera_map_text)
    model_path = None
    if model is not None:
        model_path = root / "model.json"
        common.command(
            compile_model(
                model_path, model, version, seed=compile_seed, zenith_deg=zenith_deg
            ),
            root / "compile.log",
        )
    return Chain(root, source, band, cut, camera_map, model_path)
