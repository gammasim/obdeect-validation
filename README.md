# obdeect validation

Compare each obdeect stage (stream, propagation, sources, ray tracing) and their
composition with sim_telarray, CORSIKA and LightEmission on the same physical case.
A case passes only on its scientific checks; a working harness is not a validation.

From the workspace root (checkouts `obdeect`, `simtools`, `sim_telarray`,
`simulation-models-dev` beside this repository), run one suite in the container:

```sh
bash obdeect-validation/containers/run-validation.sh full obdeect-validation/results/full
```

Suites: `smoke`, `stream`, `propagation`, `sources`, `optics`, `chain`, `curvature`, `endtoend`, `full`,
`performance` (report only). The launcher builds obdeect and simtools from the
checkouts on a digest-pinned simtools-prod image and mounts this repository
(`/validation`), sim_telarray (`/reference`) and the models (`/models`) read-only.

Shower inputs are not in Git: copy them to `inputs/` or regenerate them from
`fixtures/corsika/*.card`; `fixtures/inputs.json` lists their SHA-256.

For harness development (Python >=3.14):

```sh
python -m pip install -e '.[dev,references]'
python -m pytest
obdeect-validate list --suite suites/smoke.json --json
obdeect-validate run --case cases/chain/optics-bypass1.json --environment environments/workspace.json --output results/one --non-interactive
```

Former `obdeect/validation/` scripts, now here: the flasher, LED and stream cross-checks are
cases (`cases/sources/flasher-dish.json`, `flasher-secondary.json`, `led-reference.json`,
`cases/stream/corsika-roundtrip.json`). The obdeect chain test runs without the container
once obdeect is installed (`python -m pytest tests/test_obdeect_chain.py`; it skips otherwise).
Summary of an IACT file: `python -m obdeect_validation.adapters.summarize_iact FILE.iact.gz`.

Case files are JSON. Each has an id, executable requirements, commands (argument
arrays, without shell expansion), and metrics with metric/expected/atol.
Commands run inside an isolated case output directory and write metrics.json.
Adapters also write details.json and PNG figures; `report.html` links them per case.

Reports: `python scripts/collect_report.py RUN reports/DATE` copies verdicts, metrics,
details and figures of a run and writes `reports/DATE/README.md`, the page to read: overview
table, then one section per validation step (stream, propagation, sources, optics, chain)
with each case's verdict, headline numbers and captioned figures. Performance results come
last, as their own section with an overview table and per-stage tables. Titles,
headline numbers and captions are set in `scripts/report_layout.json`. Hand-written text
(run, versions, findings) goes into the `<!-- notes:KEY -->` blocks of the README; they are
kept when it is rewritten (`--readme-only reports/DATE` rewrites without copying).
Collecting replaces the report (case directories of an earlier collection are removed; the
README notes stay). A separate run, such as performance, is added with
`--add RUN reports/DATE`: its cases are listed in `report_only.json` and shown as REPORT ONLY.
Figures use one palette (`adapters/plotting.py`): obdeect blue, reference program
orange, analytic expectation black dashed, one colour per telescope or zenith angle.
Shared helpers (commands, sim_telarray, photon readers, output files) are in
`adapters/common.py` and `adapters/readers.py`. `VALIDATION_MODELS` overrides the
simulation-models location (default `/models`).
A case may list `inputs` (absolute file paths, for example the shower files under
`/validation/inputs/`): if one is missing the case is reported as blocked (exit status 3) and
its commands are not run. `inputs` of the environment file block every case.
Suites select case files. Keep inputs and seeds explicit. Never change a reference
or tolerance simply to make a failing comparison pass.

Agent loop: inspect case and reference assumptions; run doctor; run one affected
case; inspect metrics and logs; fix implementation; rerun that case and harness
tests. Report commands and results, including blocked cases. Do not alter the
reference simulator, publish results, or update scientific limits implicitly.

The implementation plan and current scientific coverage are recorded in the
workspace project-notes/projects/cta/202611/full-simulation-chain/{plan.md,pr.md}.

## Coverage

| Group | Cases | Reference |
|---|---|---|
| stream | LED to EventIO; CORSIKA round trip incl. stdin/tee variants; EventIO headers of the 16 physics-matrix shower inputs vs CORSIKA logs | eventio reader, sim_telarray camera, CORSIKA logs |
| propagation | analytic Beer–Lambert; sim_telarray tables; refraction/time | sim_telarray tables, compiled CORSIKA atmo |
| curvature | Chapman extinction; spherical geometry; curved Bouguer refraction; input guards (no showers); curved flag, emission geometry and curved extinction of CORSIKA CURVED showers; high-zenith chain | analytic, independent ray tracer, CORSIKA CURVED, sim_telarray (plane-parallel) |
| sources | LED, flasher dish/secondary, point source, beam scan, laser, night sky, muon; illuminators ILLN-01/ILLS-01 at 3 distances; night-sky camera p.e.; scattered laser beam | pixled, ff-1m, ff-gct, xyzls, xyls, ls-beam, nsbls, fake-muon, analytic expectation (+ camera) |
| optics | LST, MST, SST: D80, effective area, focal length, 0–4 deg | simtools-validate-optics, both backends |
| chain | near-threshold round trip; obdeect atmosphere; obdeect optics + Bypass_Optics=1 (LST, MST, SST; 40/60 deg, azimuth 180, proton, South); mode-3 envelope (Bypass_Optics=3, QE envelope) | direct sim_telarray |
| endtoend | photon-identical output across threads, chunking and pipes; tee fan-out to two models; tampered efficiency ledgers rejected; damaged streams fail cleanly | exact equality, separate full runs, analytic controls |
| performance | air-shower array and proton, pipe timing (report only) | sim_telarray on the same input |

Known FAIL: `sources.illuminator-illn-01` (`half_cx_mean` limit, one of three distances). It fails the same way on the main obdeect build, so it is not caused by the stream or propagation changes; the cause is open.

### Curvature suite

`suites/curvature.json` validates the curved and spherical atmosphere of the obdeect
propagation stage. Four cases need no showers (`propagation.curvature-chapman`,
`propagation.curvature-geometry`, `propagation.refraction-curved-ray`,
`propagation.curvature-input-guards`) and use analytic or independent references. Four need
CORSIKA CURVED showers at 65 and 75 deg (`stream.corsika-curved-flag`,
`propagation.curvature-corsika-geometry`, `propagation.curvature-corsika-extinction`,
`chain.curvature-high-zenith`); the cards are in `fixtures/corsika/*-curved.card` and the
files are generated on a host with the curved CORSIKA binary (variants in
`$SIMTOOLS_CORSIKA_PATH/build_opts.yml`) with
`python -m obdeect_validation.adapters.curved_showers run --geometry curved --zenith 65 ...`,
then copied to `inputs/` (hashes in `fixtures/inputs.json`). The run
checks the CURVED banner of the CORSIKA log, so a flat binary cannot pass as curved.
sim_telarray applies a plane-parallel airmass 1/|cz| also to CURVED input, so the chain case
tests the predicted curved-to-plane-parallel ratio, not a zero difference. The comparison with the
CORSIKA `DISTIP`/`TOFIP` tables needs an instrumented CORSIKA run, see
`docs/curvature-corsika-tables.md`. Details of all curvature cases: `docs/curvature-validation.md`.

Optics limits (2% D80, 1% area, 0.2% focal length) are approved 2026-11-15, taken from the
largest deviations measured in the old validation.
Current verdicts and open issues are in the workspace project notes.

Exit status: 0 pass, 1 scientific difference, 2 execution/configuration error,
3 missing prerequisite. A blocked test does not count as a pass.
