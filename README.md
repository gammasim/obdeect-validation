# obdeect validation

Compare each obdeect stage (stream, propagation, sources, ray tracing) and their
composition with sim_telarray, CORSIKA and LightEmission on the same physical case.
A case passes only on its scientific checks; a working harness is not a validation.

From the workspace root (checkouts `obdeect`, `simtools`, `sim_telarray`,
`simulation-models-dev` beside this repository), run one suite in the container:

```sh
bash obdeect-validation/containers/run-validation.sh full obdeect-validation/results/full
```

Suites: `smoke`, `stream`, `propagation`, `sources`, `optics`, `optics-window`, `chain`, `curvature`, `endtoend`, `full`,
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

The obdeect-only regression suite runs on a laptop without containers, through the same runner
and the host environment file `environments/host.json`. The file sets no paths: it requires `obdeect`,
`python3` and `cc` on `PATH` and the Python modules the cases import, and the runner inherits the
caller's environment. Export `SIMTEL_REFERENCE_SOURCE`, `VALIDATION_ROOT` and `VALIDATION_MODELS`
first where a case needs them (sim_telarray and simulation-models checkouts); `refraction-high-zenith`
also needs the CORSIKA atmosphere and hessio sources (`CORSIKA_ATMO_DIR`, `CORSIKA_SUPPORT_DIR`,
`CORSIKA_INCLUDE_DIRS`, see its adapter). Then:

```sh
obdeect-validate run --suite suites/regression.json --environment environments/host.json \
  --output results/regression --non-interactive
```

Former `obdeect/validation/` scripts, now here: the flasher, LED and stream cross-checks are
cases (`cases/sources/flasher-dish.json`, `flasher-secondary.json`, `led-reference.json`,
`cases/stream/corsika-roundtrip.json`). The obdeect chain test runs without the container
once obdeect is installed (`python -m pytest tests/test_obdeect_chain.py`; it skips otherwise).
Summary of an IACT file: `python -m obdeect_validation.adapters.summarize_iact FILE.iact.gz`.

The LST camera-window study is `suites/optics-window.json`. It is a study, not a gated validation:
its three cases carry no metrics (their names say so) and the suite is not part of `suites/full.json`.
It needs obdeect with camera-window tracing (obdeect branch `lst-window`, not yet in obdeect `main`)
and the LSTN-design window parameter records (`camera_window_parameters`,
`camera_window_refractive_index`, `camera_window_internal_transmittance`, version 1.0.0), which are
not in the simulation-models production checkout; the study stops with an error while they are
missing. Its config sets `simulation_models_path` to `/models` (container). Run it with:

```sh
obdeect-validate run --suite suites/optics-window.json \
  --output results/window-study --environment environments/workspace.json --non-interactive
```

It runs obdeect only:
window off at 6.55 cm, window on at 6.55 cm, then window on from 6.95 to 6.15 cm
in 0.05 cm steps. The scan reports its best on-axis D80 and compares its 1-4 deg
D80 values to the window-off case. The study adds the window records to a temporary copy of the
production table, so the checkout is not changed.

Case files are JSON. Each has an id, executable requirements, commands (argument
arrays, without shell expansion), and metrics with metric/expected/atol. Case `argv` and `inputs`
may use the placeholders `{python}` (this interpreter), `{validation}` (`VALIDATION_ROOT`, default
`/validation`) and `{reference}` (`SIMTEL_REFERENCE_SOURCE`, default `/reference`). The `env` block
of the environment file overrides the process environment; the container mounts keep the defaults,
so the host sets `VALIDATION_ROOT` and `SIMTEL_REFERENCE_SOURCE` before a run (see `environments/host.json`).
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

## Performance run (timing, report only)

Needs: Linux or macOS with podman; the workspace checkouts `obdeect`, `simtools`, `sim_telarray`
(branch `bypass-obdeect-mode-3`: the candidate chain needs sim_telarray `Bypass_Optics=3`),
`simulation-models-dev`, and a **git clone of hessio** that contains commit
`604578331f4696703001e4af3b95a19cdcaa3852` (`HESSIO_REPO`; a CTAO GitLab clone needs credentials).
The images are not in a registry: the first call builds them from the checkouts (`BUILD_ONLY=1`),
later calls reuse them with `SKIP_BUILD=1`.

Base image: `run-validation.sh` pins `simtools-prod:v0.39.0-v78010-...@sha256:...` (CORSIKA 7.8010,
`corsika_qgs2_urqmd_flat`); a newer local simtools-prod image may lack that CORSIKA. The CORSIKA
interaction tables are not in the image: clone `corsika7-interaction-tables` and use its
`interaction-tables` directory as `CORSIKA_DATA` below.

Shower inputs (not in Git; `fixtures/inputs.json` lists their SHA-256, about 2.6 GB for the two
100-shower files): `inputs/gamma-20deg-100-3000GeV-100showers.iact`,
`inputs/gamma-20deg-100-3000GeV-array-north-100showers.iact`, `inputs/proton-20deg-300-10000GeV.iact`
(45 showers), each with its `.corsika.log`. Generate a missing one from its card in `fixtures/corsika/`:

```sh
CORSIKA_DATA=/path/to/corsika7-interaction-tables/interaction-tables \
  bash obdeect-validation/containers/generate-input.sh gamma-20deg-100-3000GeV-100showers   # a few minutes
```

CORSIKA reproduces the showers (same bunches and energies on other machines), but the file header
carries a date, so the SHA-256 of a regenerated file differs; compare the `.corsika.log` instead.
The hashes in `fixtures/inputs.json` are those of the files in `inputs/` regenerated on 2026-11-23 and
2026-11-24 (per-file `generated` date); the earlier hashes are kept as `previous_sha256`.
`generate-input.sh` mounts sim_telarray at `/reference` and the CORSIKA interaction tables at
`/corsika-data`; CORSIKA reads its tables from there (`DATDIR /corsika-data/`), so the tables keep
their original names. (`generate-input.sh` has not been run on Linux.)

Disk: images about 2 GB each (two), inputs 2.6 GB, each case a few hundred MB while it runs, the
results 20 MB per run. On macOS the podman VM does not give freed space back to the host; run
`podman machine ssh "sudo fstrim -a"` after a run.

Timing needs a quiet machine: stop other containers and check `uptime` (load below 1). Run the
cases one after another, each on its own, pinned to one core (`CPUS=0`), so that wall and CPU
time are comparable and the pipe of five programs shares one core:

```sh
export MODE3=1 ENVIRONMENT=performance HESSIO_REPO=/path/to/hessio CPUS=0
out=obdeect-validation/results/perf-$(date +%F)
BUILD_ONLY=1 bash obdeect-validation/containers/run-validation.sh performance $out   # builds the images, runs nothing
export SKIP_BUILD=1
for c in air-shower air-shower-proton air-shower-array bottleneck stages scaling stream-pipe flasher-dish point-source; do
  CASE=cases/performance/$c.json bash obdeect-validation/containers/run-validation.sh performance $out/$c
done
```

The nine cases together take about 30 minutes on one core, the array case about 15 of them. Without
`CASE` a call runs all nine cases.

`CPUS=0,1` would give two cores. Without `CPUS` the container sees every core of the host (or of
the podman VM on macOS). To collect the report: `python scripts/collect_report.py --add $out/<case> reports/DATE`
for each case (adds them as REPORT ONLY), then fill the `notes:performance*` blocks of
`reports/DATE/README.md`. The container's `container.json` records the images, the sim_telarray
commit and the host; the `host` line of the report comes from it.

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

Open validation design issue: `chain.optics-bypass1-40deg` (100-3000 GeV) and `chain.optics-bypass1-60deg`
(300-5000 GeV) use different energy ranges, so their zenith comparison is confounded by energy and the
zenith dependence is not isolated. The physics inputs are unchanged (the names say so); a matched-energy
pair needs new shower inputs.

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

Optics requirement (set by the project, 2026-11-24): D80 and effective area agree with sim_telarray
to 1 % (goal) and at most 2 % (acceptance, gated) per offset; effective focal length 0.2 %. The
limits are not widened by statistical errors: a comparison whose 3-sigma error exceeds the limit is
reported as underpowered and fails.
Current verdicts and open issues are in the workspace project notes.

Exit status: 0 pass, 1 scientific difference, 2 execution/configuration error,
3 missing prerequisite. A blocked test does not count as a pass.
