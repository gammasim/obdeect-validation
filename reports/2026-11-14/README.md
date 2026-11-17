# obdeect versus sim_telarray, 2026-11-14

<!-- notes:intro -->
Run: `bash obdeect-validation/containers/run-validation.sh full obdeect-validation/results/full-7.0.0`
(suite `full`, 19 cases) in the simtools-prod based container; obdeect `0.8.1.dev52+g6721ba04b`,
sim_telarray v2025-11-30-rc, simulation models 7.0.0 for all cases (optics and chain). Wall time
908 s. Supersedes the earlier 2026-11-14 run that used models 6.0.2 for optics.

Change against the 6.0.2 run (same cases, same limits):

- Optics: numbers differ per telescope; all within the provisional limits (LST D80 0.53 %, MST
  1.32 %, SST area 0.33 % vs 0.39 % before).
- Stream, chain (atmosphere and optics) and propagation: unchanged within the printed precision,
  except the sim_telarray table error (0 → 4.5e-8, still well inside tolerance).
- Sources: camera p.e. ratios moved by about 0.5–1 % (flasher dish 1.0019 → 1.0013; point source
  0.999 → 1.010; flasher secondary 0.996 → 1.001; muon track 0.981 → 0.997). Night sky unchanged
  (0.80). The cause was not isolated: obdeect moved from `0.8.1.dev31` to `dev52` as well as the model.
<!-- /notes -->

**Result: 19 PASS** (suite `full`, status PASS). Each case below shows what is compared, its headline numbers and figures; all metrics are in the case's `metrics.json`, limits and settings in `details.json`.

[6. Performance](#step-performance) is at the bottom: 3 FAIL, 6 PASS.

| Step | Case | Headline | Verdict |
|---|---|---|---|
| 1. Photon stream format | [LED stream → sim_telarray EventIO, strict reader](#stream-led-eventio) | photons 60 / 60, block order as CORSIKA: 1 | PASS |
| 1. Photon stream format | [CORSIKA → stream → sim_telarray](#stream-corsika-roundtrip) | p.e. +0.07% ± 0.09%, 0 trigger mismatches | PASS |
| 2. Atmosphere (propagation) | [Extinction vs analytic Beer–Lambert](#propagation-extinction-analytic) | largest transmission error 0.0e+00 | PASS |
| 2. Atmosphere (propagation) | [Extinction vs sim_telarray tables](#propagation-extinction-simtel-table) | largest error 4.5e-08 | PASS |
| 2. Atmosphere (propagation) | [Refraction and timing vs CORSIKA atmo](#propagation-atmosphere-corsika) | displacement 4.4e-04 m, time 1.6e-03 ns | PASS |
| 3. Light sources | [LED vs pixled](#sources-led-reference) | photons 60 / 60 | PASS |
| 3. Light sources | [Flasher on the dish vs ff-1m](#sources-flasher-dish) | photons 1.0015 ± 0.0008, camera p.e. 1.0013 ± 0.0019 | PASS |
| 3. Light sources | [Flasher via secondary mirror vs ff-gct](#sources-flasher-secondary) | photons 1.0002 ± 0.0022, time +0.008 ns, camera p.e. 1.001 ± 0.006 | PASS |
| 3. Light sources | [Point source vs xyzls](#sources-point-source) | accepted photons 0.995 ± 0.012, time +0.004 ns, camera p.e. 1.010 ± 0.010 | PASS |
| 3. Light sources | [Beam scan vs xyls](#sources-beam-scan) | photons 1.000 ± 0.006 | PASS |
| 3. Light sources | [Laser beam geometry vs ls-beam](#sources-laser-beam) | path 131.60 vs 131.65 m | PASS |
| 3. Light sources | [Night sky vs analytic and nsbls](#sources-night-sky) | rate error 9.2e-08, area-uniform KS 0.0006, camera p.e. vs nsbls 0.80 (recorded only) | PASS |
| 3. Light sources | [Muon track vs fake-muon and analytic](#sources-muon-track) | yield 1.0003 ± 0.0029, time vs analytic -1.1e-05 ns, camera p.e. 0.997 ± 0.009 | PASS |
| 4. Telescope optics | [Optics lst-7.0.0](#optics-lst-7-0-0) | D80 ≤ 0.53%, area ≤ 0.09%, focal length ≤ 0.018% | PASS |
| 4. Telescope optics | [Optics mst-7.0.0](#optics-mst-7-0-0) | D80 ≤ 1.32%, area ≤ 0.18%, focal length ≤ 0.059% | PASS |
| 4. Telescope optics | [Optics sst-7.0.0](#optics-sst-7-0-0) | D80 ≤ 0.13%, area ≤ 0.33%, focal length ≤ 0.010% | PASS |
| 5. Full chain | [Near-threshold showers (20–100 GeV), stream round trip](#chain-stream-roundtrip-threshold) | 291 telescope events, p.e. +0.00%, 0 trigger mismatches | PASS |
| 5. Full chain | [obdeect atmosphere + sim_telarray optics and camera](#chain-atmosphere-substitution) | p.e. +0.12% ± 0.15%, image χ² 1.12 | PASS |
| 5. Full chain | [obdeect atmosphere and optics + sim_telarray camera](#chain-optics-bypass1) | p.e. +0.06% ± 0.19%, image χ² 0.99, 0 trigger mismatches | PASS |

## How to read this report

**What is compared.** Each case compares obdeect with a reference: an analytic expectation, a sim_telarray, LightEmission or CORSIKA program, or both. Signed values are obdeect / reference − 1, ratios are obdeect / reference, ± is one standard error.

**What PASS means.** All gated checks of the case are inside limits fixed before the run (`details.json`). Where a reference program is known to be wrong, the case is gated against the analytic expectation, and the comparison with that program is *recorded only*: shown, but not part of the verdict.

**Same frame before comparing.** obdeect light sources write photons where they are emitted, the LightEmission programs where they arrive. Before comparing, obdeect photons are moved to the reference plane (flight time added). The time panels show the arrival times as compared.

**Known features of the reference programs.** Some differences in the plots come from the reference, not from obdeect; each affected case has an *Expected differences* box under its figure.

- ff-1m, ff-gct and fake-muon write times in steps of 1/32 ns: their time histograms show spikes or a comb, obdeect times are continuous. Means and widths are compared.
- nsbls samples the photon entry point wrongly (radius `1 − sqrt(u)` instead of `sqrt(u)`, rotated azimuth): its entry-point, time and camera distributions differ from the analytic expectation, which obdeect follows.
- fake-muon times deviate slightly from the analytic distribution; obdeect follows it.

**Image χ²/ndf (chain cases).** Camera images of the same event and telescope from the two runs, pixel by pixel: mean over pixels with signal of (a − b)² / (a + b), with a, b the photoelectrons of the pixel in each run. Differences from Poisson fluctuations alone give about 1; systematic shape differences give more (limit 2). Values below 1 mean the two runs share part of their fluctuations (e.g. identical photons in the stream round trip). For light sources the summed images are compared the same way (χ² per pixel of the normalised images).

**KS.** Largest difference of two cumulative distributions, compared with the critical value for the sample sizes (in brackets where shown).

**Optics.** The star is moved off axis along x (E/W) and along y (N/S). The optics are not rotationally symmetric, so off axis the two directions give different values; each figure shows both, and both programs must agree for each direction.

## 1. Photon stream format

The obdeect photon stream must carry CORSIKA photons and light-source photons to sim_telarray without loss: CORSIKA → stream → sim_telarray EventIO, compared with sim_telarray on the original CORSIKA file.

<!-- notes:stream -->

<!-- /notes -->

<a id="stream-led-eventio"></a>

### LED stream → sim_telarray EventIO, strict reader — PASS

Case: LED light-source stream to sim_telarray EventIO, read by the independent eventio reader

An LED light-source stream is converted to sim_telarray EventIO and read back with the independent `eventio` reader.

| Quantity | Value |
|---|---|
| Events (stream / EventIO) | 3 / 3 |
| Photons (stream / EventIO) | 60 / 60 |
| Block order as CORSIKA (1 = yes) | 1 |

`stream.led-eventio`: 9/9 gated checks pass. [metrics.json](stream.led-eventio/metrics.json) · [details.json](stream.led-eventio/details.json)

<a id="stream-corsika-roundtrip"></a>

### CORSIKA → stream → sim_telarray — PASS

Case: CORSIKA IACT -> obdeect stream -> sim_telarray IACT, bunch and camera comparison, plus stream command variants

CORSIKA IACT file → `obdeect stream from-corsika` → `obdeect stream to-simtel` → sim_telarray, against sim_telarray on the CORSIKA file. Also checks the stdin, cat and tee variants. Array: three PROD6 LSTs (sim_telarray `CTA-PROD6-LST.cfg`, not the simulation-models release) at the CTAO-North altitude of 2156 m, at x, y = (0, 0), (50, 0) and (0, 50) m; gamma showers at 20° zenith.

| Quantity | Value |
|---|---|
| Telescope events compared | 300 |
| Camera p.e. bias | +0.07% ± 0.09% |
| Trigger mismatches | 0 |
| Bunch count difference | 0 |
| Largest codec error (in units of the codec bound) | 0.076 |

`stream.corsika-roundtrip`: 14/14 gated checks pass. [metrics.json](stream.corsika-roundtrip/metrics.json) · [details.json](stream.corsika-roundtrip/details.json)

![Camera p.e. per telescope event, colour per telescope (type and position in the legend): round trip against direct sim_telarray (left), relative difference (right).](stream.corsika-roundtrip/camera-comparison.png)

*Camera p.e. per telescope event, colour per telescope (type and position in the legend): round trip against direct sim_telarray (left), relative difference (right).*

## 2. Atmosphere (propagation)

obdeect propagation (extinction, refraction, travel time) against analytic results, sim_telarray's extinction tables and the compiled CORSIKA atmosphere code.

<!-- notes:propagation -->

<!-- /notes -->

<a id="propagation-extinction-analytic"></a>

### Extinction vs analytic Beer–Lambert — PASS

Case: Beer-Lambert constant optical depth (analytic reference)

Atmospheric transmission for an exponential atmosphere against the analytic Beer–Lambert result, zenith 0–60°.

| Quantity | Value |
|---|---|
| Samples | 16 |
| Largest transmission error | 0.0e+00 |

`propagation.extinction-analytic`: 2/2 gated checks pass. [metrics.json](propagation.extinction-analytic/metrics.json)

![Transmission of obdeect and the analytic expectation.](propagation.extinction-analytic/extinction-analytic.png)

*Transmission of obdeect and the analytic expectation.*

<a id="propagation-extinction-simtel-table"></a>

### Extinction vs sim_telarray tables — PASS

Case: Real sim_telarray transmission table, telescope path and thinning

Extinction from a sim_telarray optical-depth table against sim_telarray's own interpolation, for wavelengths, heights, zenith angles, telescope elevations and in-telescope path lengths; weighted and thinned modes.

| Quantity | Value |
|---|---|
| Grid (wavelengths × heights × zeniths) | 9 × 5 × 4 |
| Largest error (weighted mode) | 4.5e-08 |
| Largest error (in-telescope depth) | 4.3e-19 |
| Thinned count residual | 2.77 σ |

`propagation.extinction-simtel-table`: 10/10 gated checks pass. [metrics.json](propagation.extinction-simtel-table/metrics.json)

![Weighted transmission of obdeect against the sim_telarray table.](propagation.extinction-simtel-table/extinction-weighted.png)

*Weighted transmission of obdeect against the sim_telarray table.*

<a id="propagation-atmosphere-corsika"></a>

### Refraction and timing vs CORSIKA atmo — PASS

Case: Atmosphere refraction against compiled CORSIKA IACT/atmo

Refraction and travel time against the compiled CORSIKA atmosphere code (atmo.c) for emission heights 3–20 km and zenith 0–60°.

| Quantity | Value |
|---|---|
| Samples (bent) | 24 (24) |
| Largest displacement error | 4.4e-04 m |
| Largest travel-time error | 1.6e-03 ns |
| Largest direction-cosine error | 1.6e-12 |

`propagation.atmosphere-corsika`: 5/5 gated checks pass. [metrics.json](propagation.atmosphere-corsika/metrics.json)

![obdeect − CORSIKA: direction cosine, displacement and travel time against emission height, per zenith angle.](propagation.atmosphere-corsika/atmosphere-reference.png)

*obdeect − CORSIKA: direction cosine, displacement and travel time against emission height, per zenith angle.*

## 3. Light sources

obdeect light sources against the sim_telarray LightEmission programs and analytic expectations: photon yield, directions, times and wavelengths at emission, then camera photoelectrons in sim_telarray. Signed values are obdeect / reference − 1; ratios are obdeect / reference.

<!-- notes:sources -->
Reference findings (recorded, not judged; details in each `details.json`):

- fake-muon: arrival-time distribution differs from the analytic one by KS 0.021 (mean −5 ps).
- nsbls: entry radius `1 − sqrt(u)` instead of `sqrt(u)` (flux not isotropic, too much light
  at the dish centre), and the entry point is rotated with the azimuth of `d` instead of `−d`
  (`IactLightEmission.cc` ~l.1395), which breaks impact parameters and times. The night-sky camera
  ratio (0.80) is therefore recorded only, not gated; with the same sampler in both programs
  (2026-11-13 run) it was 0.998 ± 0.009.
- ff-gct: near-horizontal photons (|cz| < 0.3) carry a refraction time term of up to 4 ns.
- sim_telarray skips an event block larger than its input buffer (default limits) and still
  exits 0; the flasher-secondary camera stage therefore uses 1e6 photons per flasher. Not an
  obdeect format problem: with `CORSIKA_IO_BUFFER=2GiB MAX_BUNCHES=50000000` the 1.28 GB
  stream is read completely (568833 p.e., 10x the 1e6 run).
- Point source: the raw photon ratio is recorded but not gated; the gated quantity is the
  accepted-photon ratio.
<!-- /notes -->

<a id="sources-led-reference"></a>

### LED vs pixled — PASS

Case: sources.led-reference

LED in the camera against the pixled reference program: identical photons expected.

| Quantity | Value |
|---|---|
| Photons (obdeect / reference) | 60 / 60 |
| Largest per-event photon difference | 0.0 |
| Position centroid difference | 1.05e-08 cm |

`sources.led-reference`: 8/8 gated checks pass. [metrics.json](sources.led-reference/metrics.json) · [details.json](sources.led-reference/details.json)

<a id="sources-flasher-dish"></a>

### Flasher on the dish vs ff-1m — PASS

Case: Dish flasher versus ff-1m (emission; camera stage in details.json)

Flasher on the dish axis against ff-1m: emission (stage 1), then camera photoelectrons through sim_telarray (stage 2).

| Quantity | Value |
|---|---|
| Photons per event | 1.0015 ± 0.0008 |
| Direction mean cx / cy | +0.62 σ / -0.66 σ |
| Arrival time mean (after frame correction) | -0.006 ± 0.000 ns |
| Camera p.e. | 1.0013 ± 0.0019 |
| Camera image χ²/ndf | 1.03 |

`sources.flasher-dish`: 14/14 gated checks pass. [metrics.json](sources.flasher-dish/metrics.json) · [details.json](sources.flasher-dish/details.json)

![Emission: direction cosines, arrival time at the reference plane (obdeect times moved to the plane, as in the metrics) and wavelength of obdeect and LightEmission photons.](sources.flasher-dish/stage1.png)

*Emission: direction cosines, arrival time at the reference plane (obdeect times moved to the plane, as in the metrics) and wavelength of obdeect and LightEmission photons.*

> **Expected differences.** Time panel: ff-1m writes times in steps of 1/32 ns, so its histogram is a row of spikes; obdeect times are continuous. Both cover the same window; the mean difference is -0.006 ns (limit 0.05 ns). Wavelength: both are monochromatic.

<a id="sources-flasher-secondary"></a>

### Flasher via secondary mirror vs ff-gct — PASS

Case: Secondary flasher versus ff-gct (stage 1 emission)

Flashers at the camera corners, light reflected by the secondary mirror, against ff-gct; camera photoelectrons through sim_telarray.

| Quantity | Value |
|---|---|
| Photons per event | 1.0002 ± 0.0022 |
| Arrival time mean | +0.008 ± 0.003 ns |
| Camera p.e. | 1.0014 ± 0.0059 |
| Camera image χ²/ndf | 1.08 |

`sources.flasher-secondary`: 12/12 gated checks pass. [metrics.json](sources.flasher-secondary/metrics.json) · [details.json](sources.flasher-secondary/details.json)

![Emission: direction cosines, arrival time at the reference plane (obdeect times moved to the plane, as in the metrics) and wavelength of obdeect and LightEmission photons.](sources.flasher-secondary/stage1.png)

*Emission: direction cosines, arrival time at the reference plane (obdeect times moved to the plane, as in the metrics) and wavelength of obdeect and LightEmission photons.*

> **Expected differences.** Time panel: ff-gct writes times in steps of 1/32 ns, so its histogram has regular spikes on the same falling distribution; obdeect times are continuous. Mean difference +0.008 ns (limit 0.05 ns). Times are compared with the ff-gct time model: photons with |cz| ≥ 0.3 within the camera radius on both sides.

<a id="sources-point-source"></a>

### Point source vs xyzls — PASS

Case: Point source versus xyzls (emission; camera stage in details.json)

Isotropic point source at a distance against xyzls; photons accepted by the telescope and camera photoelectrons.

| Quantity | Value |
|---|---|
| Accepted photons (gated) | 0.995 ± 0.012 |
| Raw photons (recorded only) | 1.039 ± 0.001 |
| Arrival time mean | +0.004 ± 0.001 ns |
| Camera p.e. | 1.010 ± 0.010 |

`sources.point-source`: 14/14 gated checks pass. [metrics.json](sources.point-source/metrics.json) · [details.json](sources.point-source/details.json)

![Emission: direction cosines, arrival time at the reference plane (obdeect times moved to the plane, as in the metrics) and wavelength of obdeect and LightEmission photons.](sources.point-source/stage1.png)

*Emission: direction cosines, arrival time at the reference plane (obdeect times moved to the plane, as in the metrics) and wavelength of obdeect and LightEmission photons.*

<a id="sources-beam-scan"></a>

### Beam scan vs xyls — PASS

Case: Beam scan versus xyls (emission; camera stage in details.json)

Scanning laser beam against the LightEmission beam program.

| Quantity | Value |
|---|---|
| Photons per event | 1.000 ± 0.006 |
| Direction mean cx / cy | +0.00 σ / +0.00 σ |

`sources.beam-scan`: 12/12 gated checks pass. [metrics.json](sources.beam-scan/metrics.json) · [details.json](sources.beam-scan/details.json)

![Emission: direction cosines, arrival time at the reference plane (obdeect times moved to the plane, as in the metrics) and wavelength of obdeect and LightEmission photons.](sources.beam-scan/stage1.png)

*Emission: direction cosines, arrival time at the reference plane (obdeect times moved to the plane, as in the metrics) and wavelength of obdeect and LightEmission photons.*

<a id="sources-laser-beam"></a>

### Laser beam geometry vs ls-beam — PASS

Case: Laser beam versus ls-beam: FoV crossing geometry only (yield, direction, time, camera BLOCKED)

Laser beam geometry against ls-beam: where the beam enters and leaves the telescope field of view.

| Quantity | Value |
|---|---|
| Path in the field of view (obdeect / reference) | 131.60 / 131.65 m |
| Entry / exit offset | 0.00 / 0.05 m |

`sources.laser-beam`: 6/6 gated checks pass. [metrics.json](sources.laser-beam/metrics.json) · [details.json](sources.laser-beam/details.json)

![Beam path of obdeect and ls-beam relative to the field of view.](sources.laser-beam/laser-beam.png)

*Beam path of obdeect and ls-beam relative to the field of view.*

<a id="sources-night-sky"></a>

### Night sky vs analytic and nsbls — PASS

Case: Night sky versus nsbls: yield, rate, direction, projected area, time, camera p.e.

Night-sky background source against the analytic flux and nsbls; camera photoelectrons through sim_telarray.

| Quantity | Value |
|---|---|
| Rate, relative error vs configured | 9.2e-08 |
| Area-uniform entry (KS, critical) | 0.0006 (0.0011) |
| Flat pulse (KS, critical) | 0.0006 (0.0011) |
| Camera p.e. vs nsbls (recorded only) | 0.796 ± 0.007 |

`sources.night-sky`: 12/12 gated checks pass. [metrics.json](sources.night-sky/metrics.json) · [details.json](sources.night-sky/details.json)

![Angular distribution in the cone, projected entry point (uniform expected), time at the centre plane, and camera p.e. per triggered event.](sources.night-sky/night-sky.png)

*Angular distribution in the cone, projected entry point (uniform expected), time at the centre plane, and camera p.e. per triggered event.*

> **Expected differences.** The nsbls curves differ because nsbls is wrong, not obdeect. Projected entry point (top right): an isotropic flux is flat in ρ²; obdeect is flat (KS 0.0006, critical 0.0011), nsbls piles up at the centre (KS 0.45) because it draws the radius as 1 − sqrt(u) instead of sqrt(u) and rotates the entry point with the wrong azimuth (`IactLightEmission.cc` ~l.1395). Time (bottom left): the wrong entry points give nsbls a spread of 3.6 ns; obdeect has the configured pulse (0.29 ns). Camera (bottom right): nsbls puts more light on the dish centre, so its camera sees more; obdeect / nsbls = 0.80, recorded only. With the nsbls sampler in obdeect (2026-11-13 report) the camera ratio was 0.998 ± 0.009.

<a id="sources-muon-track"></a>

### Muon track vs fake-muon and analytic — PASS

Case: Muon track versus fake-muon: yield, Cherenkov angle, spectrum, time, camera p.e.

Muon track against fake-muon and the analytic Cherenkov expectation; camera photoelectrons through sim_telarray.

| Quantity | Value |
|---|---|
| Photon yield (obdeect / fake-muon) | 1.0003 ± 0.0029 |
| Time vs analytic | -1.1e-05 ns, KS 0.0002 |
| Camera p.e. | 0.997 ± 0.009 |

`sources.muon-track`: 18/18 gated checks pass. [metrics.json](sources.muon-track/metrics.json) · [details.json](sources.muon-track/details.json)

![Wavelength against 1/λ², Cherenkov angle residual, arrival time against the analytic distribution (black), and camera p.e. per triggered event.](sources.muon-track/muon-track.png)

*Wavelength against 1/λ², Cherenkov angle residual, arrival time against the analytic distribution (black), and camera p.e. per triggered event.*

> **Expected differences.** Time (bottom left): fake-muon writes times in steps of 1/32 ns (comb); obdeect follows the analytic distribution (black; mean difference -1.1e-05 ns, KS 0.0002, critical 0.0056). Cherenkov angle (top right): the edges differ by about 1e-7 in cos θ, far below what an image resolves. Wavelength and camera p.e. agree.

## 4. Telescope optics

Ray tracing of obdeect against sim_telarray for the same telescope model, star light at 0–4° off axis along two directions: PSF (D80), effective mirror area and effective focal length.

<!-- notes:optics -->

<!-- /notes -->

<a id="optics-lst-7-0-0"></a>

### Optics lst-7.0.0 — PASS

Case: Optics LSTN-01 (model 7.0.0): simtools sim_telarray vs obdeect, offsets 0-4 deg, N/E

sim_telarray ray tracing (via simtools) against obdeect for the same telescope model: a star at 10 km, 0–4° off axis, moved along x (E/W) and along y (N/S). PSF (D80), effective mirror area and effective focal length per offset.

| Quantity | Value |
|---|---|
| Offsets | 9 |
| Largest PSF (D80) difference | 0.53% |
| Largest effective area difference | 0.09% |
| Largest effective focal length difference | 0.018% |

`optics.lst-7.0.0`: 5/5 gated checks pass. [metrics.json](optics.lst-7.0.0/metrics.json) · [details.json](optics.lst-7.0.0/details.json)

![obdeect (blue) and sim_telarray (orange) against off-axis angle; filled circles and solid lines: star offset along x (E/W), open squares and dashed lines: along y (N/S). On axis there is one point.](optics.lst-7.0.0/optics_values.png)

*obdeect (blue) and sim_telarray (orange) against off-axis angle; filled circles and solid lines: star offset along x (E/W), open squares and dashed lines: along y (N/S). On axis there is one point.*

![obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = provisional limit.](optics.lst-7.0.0/optics_relative_difference.png)

*obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = provisional limit.*

Plots written by simtools (one program each): [obdeect/validate_optics_LSTN-01_eff_area.png](optics.lst-7.0.0/obdeect/validate_optics_LSTN-01_eff_area.png) · [obdeect/validate_optics_LSTN-01_eff_flen.png](optics.lst-7.0.0/obdeect/validate_optics_LSTN-01_eff_flen.png) · [obdeect/validate_optics_LSTN-01_psf_cm.png](optics.lst-7.0.0/obdeect/validate_optics_LSTN-01_psf_cm.png) · [obdeect/validate_optics_LSTN-01_psf_deg.png](optics.lst-7.0.0/obdeect/validate_optics_LSTN-01_psf_deg.png) · [sim_telarray/validate_optics_LSTN-01_eff_area.png](optics.lst-7.0.0/sim_telarray/validate_optics_LSTN-01_eff_area.png) · [sim_telarray/validate_optics_LSTN-01_eff_flen.png](optics.lst-7.0.0/sim_telarray/validate_optics_LSTN-01_eff_flen.png) · [sim_telarray/validate_optics_LSTN-01_psf_cm.png](optics.lst-7.0.0/sim_telarray/validate_optics_LSTN-01_psf_cm.png) · [sim_telarray/validate_optics_LSTN-01_psf_deg.png](optics.lst-7.0.0/sim_telarray/validate_optics_LSTN-01_psf_deg.png)

<a id="optics-mst-7-0-0"></a>

### Optics mst-7.0.0 — PASS

Case: Optics MSTN-01 (model 7.0.0): simtools sim_telarray vs obdeect, offsets 0-4 deg, N/E

sim_telarray ray tracing (via simtools) against obdeect for the same telescope model: a star at 10 km, 0–4° off axis, moved along x (E/W) and along y (N/S). PSF (D80), effective mirror area and effective focal length per offset.

| Quantity | Value |
|---|---|
| Offsets | 9 |
| Largest PSF (D80) difference | 1.32% |
| Largest effective area difference | 0.18% |
| Largest effective focal length difference | 0.059% |

`optics.mst-7.0.0`: 5/5 gated checks pass. [metrics.json](optics.mst-7.0.0/metrics.json) · [details.json](optics.mst-7.0.0/details.json)

![obdeect (blue) and sim_telarray (orange) against off-axis angle; filled circles and solid lines: star offset along x (E/W), open squares and dashed lines: along y (N/S). On axis there is one point.](optics.mst-7.0.0/optics_values.png)

*obdeect (blue) and sim_telarray (orange) against off-axis angle; filled circles and solid lines: star offset along x (E/W), open squares and dashed lines: along y (N/S). On axis there is one point.*

![obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = provisional limit.](optics.mst-7.0.0/optics_relative_difference.png)

*obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = provisional limit.*

Plots written by simtools (one program each): [obdeect/validate_optics_MSTN-01_eff_area.png](optics.mst-7.0.0/obdeect/validate_optics_MSTN-01_eff_area.png) · [obdeect/validate_optics_MSTN-01_eff_flen.png](optics.mst-7.0.0/obdeect/validate_optics_MSTN-01_eff_flen.png) · [obdeect/validate_optics_MSTN-01_psf_cm.png](optics.mst-7.0.0/obdeect/validate_optics_MSTN-01_psf_cm.png) · [obdeect/validate_optics_MSTN-01_psf_deg.png](optics.mst-7.0.0/obdeect/validate_optics_MSTN-01_psf_deg.png) · [sim_telarray/validate_optics_MSTN-01_eff_area.png](optics.mst-7.0.0/sim_telarray/validate_optics_MSTN-01_eff_area.png) · [sim_telarray/validate_optics_MSTN-01_eff_flen.png](optics.mst-7.0.0/sim_telarray/validate_optics_MSTN-01_eff_flen.png) · [sim_telarray/validate_optics_MSTN-01_psf_cm.png](optics.mst-7.0.0/sim_telarray/validate_optics_MSTN-01_psf_cm.png) · [sim_telarray/validate_optics_MSTN-01_psf_deg.png](optics.mst-7.0.0/sim_telarray/validate_optics_MSTN-01_psf_deg.png)

<a id="optics-sst-7-0-0"></a>

### Optics sst-7.0.0 — PASS

Case: Optics SSTS-01 (dual mirror) (model 7.0.0): simtools sim_telarray vs obdeect, offsets 0-4 deg, N/E

sim_telarray ray tracing (via simtools) against obdeect for the same telescope model: a star at 10 km, 0–4° off axis, moved along x (E/W) and along y (N/S). PSF (D80), effective mirror area and effective focal length per offset.

| Quantity | Value |
|---|---|
| Offsets | 9 |
| Largest PSF (D80) difference | 0.13% |
| Largest effective area difference | 0.33% |
| Largest effective focal length difference | 0.010% |

`optics.sst-7.0.0`: 5/5 gated checks pass. [metrics.json](optics.sst-7.0.0/metrics.json) · [details.json](optics.sst-7.0.0/details.json)

![obdeect (blue) and sim_telarray (orange) against off-axis angle; filled circles and solid lines: star offset along x (E/W), open squares and dashed lines: along y (N/S). On axis there is one point.](optics.sst-7.0.0/optics_values.png)

*obdeect (blue) and sim_telarray (orange) against off-axis angle; filled circles and solid lines: star offset along x (E/W), open squares and dashed lines: along y (N/S). On axis there is one point.*

![obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = provisional limit.](optics.sst-7.0.0/optics_relative_difference.png)

*obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = provisional limit.*

Plots written by simtools (one program each): [obdeect/validate_optics_SSTS-01_eff_area.png](optics.sst-7.0.0/obdeect/validate_optics_SSTS-01_eff_area.png) · [obdeect/validate_optics_SSTS-01_eff_flen.png](optics.sst-7.0.0/obdeect/validate_optics_SSTS-01_eff_flen.png) · [obdeect/validate_optics_SSTS-01_psf_cm.png](optics.sst-7.0.0/obdeect/validate_optics_SSTS-01_psf_cm.png) · [obdeect/validate_optics_SSTS-01_psf_deg.png](optics.sst-7.0.0/obdeect/validate_optics_SSTS-01_psf_deg.png) · [sim_telarray/validate_optics_SSTS-01_eff_area.png](optics.sst-7.0.0/sim_telarray/validate_optics_SSTS-01_eff_area.png) · [sim_telarray/validate_optics_SSTS-01_eff_flen.png](optics.sst-7.0.0/sim_telarray/validate_optics_SSTS-01_eff_flen.png) · [sim_telarray/validate_optics_SSTS-01_psf_cm.png](optics.sst-7.0.0/sim_telarray/validate_optics_SSTS-01_psf_cm.png) · [sim_telarray/validate_optics_SSTS-01_psf_deg.png](optics.sst-7.0.0/sim_telarray/validate_optics_SSTS-01_psf_deg.png)

## 5. Full chain

obdeect stages inserted into the CORSIKA → sim_telarray chain, compared event by event with the direct sim_telarray run on the same showers and seeds. Bias = obdeect / reference − 1 of camera photoelectrons, with its standard error.

<!-- notes:chain -->

<!-- /notes -->

<a id="chain-stream-roundtrip-threshold"></a>

### Near-threshold showers (20–100 GeV), stream round trip — PASS

Case: Near-threshold CORSIKA showers (20-100 GeV): stream round trip and sim_telarray camera

Stream round trip for near-threshold showers (20–100 GeV), where single photons decide the trigger. Array: three PROD6 LSTs (sim_telarray `CTA-PROD6-LST.cfg`, not the simulation-models release) at the CTAO-North altitude of 2156 m, at x, y = (0, 0), (50, 0) and (0, 50) m; gamma showers at 20° zenith.

| Quantity | Value |
|---|---|
| Telescope events compared | 291 |
| Camera p.e. bias | +0.00% ± 0.00% |
| Trigger mismatches | 0 |

`chain.stream-roundtrip-threshold`: 6/6 gated checks pass. [metrics.json](chain.stream-roundtrip-threshold/metrics.json) · [details.json](chain.stream-roundtrip-threshold/details.json)

![Camera p.e. per telescope event, colour per telescope (type and position in the legend): round trip against direct sim_telarray (left), relative difference (right).](chain.stream-roundtrip-threshold/camera-comparison.png)

*Camera p.e. per telescope event, colour per telescope (type and position in the legend): round trip against direct sim_telarray (left), relative difference (right).*

<a id="chain-atmosphere-substitution"></a>

### obdeect atmosphere + sim_telarray optics and camera — PASS

Case: obdeect atmosphere substitution vs sim_telarray atmosphere (20 deg, CORSIKA band)

obdeect atmosphere (extinction) replaces sim_telarray's: CORSIKA → obdeect propagation → sim_telarray with a zero optical-depth table, against sim_telarray with its own table. Array: three PROD6 LSTs (sim_telarray `CTA-PROD6-LST.cfg`, not the simulation-models release) at the CTAO-North altitude of 2156 m, at x, y = (0, 0), (50, 0) and (0, 50) m; gamma showers at 20° zenith.

| Quantity | Value |
|---|---|
| Telescope events compared | 300 |
| Camera p.e. bias | +0.12% ± 0.15% |
| Image χ²/ndf | 1.12 |
| Trigger mismatches | 0.0 |

`chain.atmosphere-substitution`: 12/12 gated checks pass. [metrics.json](chain.atmosphere-substitution/metrics.json) · [details.json](chain.atmosphere-substitution/details.json)

![Camera p.e. per telescope event, colour per telescope (type and position in the legend): obdeect atmosphere against sim_telarray's (left), relative difference (right).](chain.atmosphere-substitution/camera-comparison.png)

*Camera p.e. per telescope event, colour per telescope (type and position in the legend): obdeect atmosphere against sim_telarray's (left), relative difference (right).*

<a id="chain-optics-bypass1"></a>

### obdeect atmosphere and optics + sim_telarray camera — PASS

Case: CORSIKA -> obdeect trace-stream -> sim_telarray camera (Bypass_Optics=1) vs direct sim_telarray, matched LSTN-01 7.0.0 model, unit atmosphere

obdeect atmosphere and optics replace sim_telarray's: CORSIKA → obdeect propagation → obdeect trace-stream → sim_telarray camera only, against the direct sim_telarray run. Telescope: one LSTN-01, simulation models 7.0.0, CTAO-North, from the configuration simtools generates; gamma showers at 20° zenith.

| Quantity | Value |
|---|---|
| Telescope events compared | 100 |
| Camera p.e. bias | +0.06% ± 0.19% |
| Image χ²/ndf | 0.99 |
| Mean arrival-time difference | +0.022 ns |
| Trigger mismatches | 0 |

`chain.optics-bypass1`: 6/6 gated checks pass. [metrics.json](chain.optics-bypass1/metrics.json) · [details.json](chain.optics-bypass1/details.json)

![Camera p.e. per telescope event, colour per telescope (type and position in the legend): obdeect chain against direct sim_telarray (left), relative difference (right).](chain.optics-bypass1/camera-comparison.png)

*Camera p.e. per telescope event, colour per telescope (type and position in the legend): obdeect chain against direct sim_telarray (left), relative difference (right).*

## Scope and limits

<!-- notes:end -->
Limits are set before the runs and recorded in each `details.json`; optics and source limits
are provisional and need scientific approval.

Not covered: `Bypass_Optics=3` (not in the sim_telarray release used here; implemented on the
sim_telarray branch `bypass-obdeect-mode-3`, in review), dual-mirror shower chain, scattered
laser light, zenith angles other than 20° in the shower chain, per-mirror and incident-angle
optics cases.
<!-- /notes -->

<a id="step-performance"></a>

## 6. Performance

Report only, apart from the timing budgets. Wall and CPU time of the obdeect chain against sim_telarray on the same input and machine, ordered from the whole chain to its parts: whole chain (gamma, proton), bottleneck inside the pipe, separate stages, fixed versus per-shower cost and threads, array, stream pipes, light sources, profiles. Timings with timing unreliable = 1 come from a loaded machine or noisy repetitions.

<!-- notes:performance -->
Run 2026-11-17 (`results/perf-2026-11-17`): image `obdeect-validation:perf` = `obdeect-validation:local` (obdeect d19279a, this validation tree) plus the mode-3 sim_telarray from `obdeect-validation:mode3`; podman VM with 6 cores, 8 GB. Other containers were running on the same VM (load average 3–7), so every timed case has timing unreliable = 1. Whole-pipe walls repeat within 1–11 %, but single stages vary by up to 83 %. Read the ratios as indicative.

- **Whole chain.** The envelope pipe takes 0.92× (gamma) and 0.90× (proton) the wall time of sim_telarray, and 1.46× and 1.38× its CPU time. The wall budgets pass and the CPU budgets fail. Physics is unchanged: the p.e. bias is within its limits in all three comparisons.
- **Where the CPU goes** in the envelope pipe (gamma, threads=1, 5.86 s against 4.00 s for the reference):

  | Stage | CPU s | Share |
  |---|---|---|
  | sim_telarray camera | 2.53 | 47 % |
  | propagation | 1.26 | 24 % |
  | trace-stream | 1.06 | 20 % |
  | from-corsika | 0.46 | 9 % |
  | to-simtel | 0.03 | 1 % |

  The camera stage costs less than in the reference (2.53 against 3.12 s), because the stream reaching it is thinned. The extra CPU comes from from-corsika and propagation, 1.72 s together, work that sim_telarray does inside its optics (0.88 s). Trace-stream on the thinned stream (1.06 s) costs about as much as the reference optics. On the unthinned stream it costs 4.1× the reference optics (gamma) and 9.8× (proton).
- **Bottleneck** (wall time each member adds to the pipe): propagation 1.26 s (37 %), sim_telarray 1.19 s (35 %), trace-stream 0.54 s (16 %), from-corsika 0.34 s (10 %), to-simtel 0.07 s (2 %). **Gain:** propagation is both the wall bottleneck and the largest CPU item without a counterpart in the reference.
- **Gains of the current design.** Thinning saves 1.08× CPU (gamma) and 1.16× (proton). Piping saves 1.46× and 1.86× against separate stages with intermediate files.
- **Fixed and per-shower cost.** The fixed costs are the model compile (0.34 s wall, 0.47 s CPU, every run, not cached) and the trace-stream start (0.24 s wall, 0.55 s CPU). Every other stage scales with the number of showers; their small negative fixed values are noise.
- **Trace-stream threads.** 2 threads give a 1.50× speedup and 4 threads 1.91×, with efficiency falling to 48 %. 8 threads were not run because the VM has 6 cores.
- **Array** (six telescopes, tee fan-out): wall 10.8 s against 45.3 s for one sim_telarray run (0.24×), CPU 24.4 s against 43.6 s (0.56×). The tee needs 0.76× the CPU of one pipe per telescope. The per-telescope p.e. bias ranges from −0.3 % to +1.1 %, all within limits. The ungated `trigger_pattern_mismatches` equals the number of matched samples: it compares the reference's six-telescope trigger pattern with a one-telescope candidate, so it means nothing here.
- **Stream pipes.** `cat | tee` itself reaches only 0.81 GB/s on this VM's file system, so the 1 GB/s budget can't be met here. `to-simtel` reaches 0.48 GB/s.
- **Light sources.** Flasher 0.70× and point source 0.78× the reference wall time; sim_telarray is almost all of it.
- **Not run.** `performance.profile`: `perf` is not in the image.
<!-- /notes -->

| Case | Headline | Verdict |
|---|---|---|
| [Air shower: obdeect chain against sim_telarray](#performance-air-shower) | envelope chain / sim_telarray: wall 0.92, CPU 1.46 (threads=1); trace-stream / sim_telarray optics CPU 4.07 | FAIL |
| [Air shower (proton): obdeect chain against sim_telarray](#performance-air-shower-proton) | envelope chain / sim_telarray: wall 0.90, CPU 1.38 (threads=1); trace-stream / sim_telarray optics CPU 9.80 | FAIL |
| [Bottleneck of the envelope pipe](#performance-bottleneck) | bottleneck propagation, pipe wall 3.4 s | PASS |
| [Stage timing (wall and CPU)](#performance-stages) | trace-stream 3.8 s, sim_telarray direct 4.3 s (wall median), timing unreliable: 1 | PASS |
| [Scaling: fixed and per-shower cost, trace-stream threads](#performance-scaling) | fixed share of trace-stream 13% | PASS |
| [Array of six telescopes: chain against one sim_telarray run](#performance-air-shower-array) | tee / per-telescope CPU 0.76, candidate tee wall 10.8 s, reference wall 45.3 s | PASS |
| [Stream pipes: throughput](#performance-stream-pipe) | to-simtel 0.48 GB/s, cat → tee 0.81 GB/s | FAIL |
| [Flasher on the dish: chain against LightEmission](#performance-flasher-dish) | candidate / reference wall 0.70 (threads=1) | PASS |
| [Point source: chain against LightEmission](#performance-point-source) | candidate / reference wall 0.78 (threads=1) | PASS |

<a id="performance-air-shower"></a>

### Air shower: obdeect chain against sim_telarray — FAIL

Case: Air shower, one LST (LSTN-01 7.0.0, 20 deg gamma): CORSIKA IACT -> sim_telarray full optics versus obdeect envelope chain (threads=1 and default) and unthinned chain; timing budgets

Reference: CORSIKA IACT → sim_telarray with full optics. Candidate: the obdeect chain up to the camera stage (optics bypassed, same photons). Envelope candidate: from-corsika, propagation (QE bound, thinned), trace-stream and to-simtel as one pipe, then sim_telarray for the camera. Run with one thread and with the obdeect default thread count. Stage shares are of the candidate's CPU time.

| Quantity | Value |
|---|---|
| Reference wall s | 4.20 |
| Reference CPU s | 4.00 |
| Reference optics CPU s | 0.88 |
| Reference camera CPU s | 3.12 |
| Candidate CPU s, threads=1 | 8.93 |
| Envelope candidate CPU s, threads=1 | 5.86 |
| Unthinned candidate CPU s, threads=1 | 6.29 |
| Trace-stream CPU / reference optics CPU | 4.07 |
| Thinning gain, CPU | 1.08 |
| Piping gain, CPU | 1.46 |
| Envelope / reference wall, threads=1 | 0.92 |
| Envelope / reference CPU, threads=1 | 1.46 |
| Timing unreliable (1 = loaded machine or noisy repetitions) | 1 |

| Stage | Candidate CPU share, threads=1 | Envelope CPU s, threads=1 | Envelope CPU share, threads=1 |
|---|---|---|---|
| from-corsika | 6% | 0.46 | 9% |
| propagation | 12% | 1.26 | 24% |
| trace-stream | 42% | 1.06 | 20% |
| to-simtel | 3% | 0.03 | 1% |
| sim_telarray | 37% | 2.53 | 47% |

| Candidate | s per shower | CPU ns per photon |
|---|---|---|
| reference | 0.21 | 215.20 |
| candidate, threads=1 | 0.44 | 479.75 |
| envelope candidate, threads=1 | 0.19 | 314.80 |

| Gated check | Actual | Limit | Verdict |
|---|---|---|---|
| `reference_completed` | 1 | 1 ± 0 | PASS |
| `candidate_threads1_completed` | 1 | 1 ± 0 | PASS |
| `candidate_threadsdefault_completed` | 1 | 1 ± 0 | PASS |
| `candidate_envelope_completed` | 1 | 1 ± 0 | PASS |
| `candidate_envelope_threads1_completed` | 1 | 1 ± 0 | PASS |
| `candidate_unthinned_threads1_completed` | 1 | 1 ± 0 | PASS |
| `physics_completed` | 1 | 1 ± 0 | PASS |
| `from_corsika_propagation_cpu_budget_exceeded` | 1 | 0 ± 0 | FAIL |
| `envelope_from_corsika_propagation_cpu_budget_exceeded` | 1 | 0 ± 0 | FAIL |
| `candidate_envelope_wall_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `candidate_envelope_cpu_budget_exceeded` | 1 | 0 ± 0 | FAIL |
| `candidate_envelope_threads1_wall_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `candidate_envelope_threads1_cpu_budget_exceeded` | 1 | 0 ± 0 | FAIL |
| `equal_statistics_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `candidate_envelope_camera_only_missing` | 0 | 0 ± 0 | PASS |
| `candidate_envelope_threads1_camera_only_missing` | 0 | 0 ± 0 | PASS |

`performance.air-shower`: 12/16 gated checks pass; failed: `from_corsika_propagation_cpu_budget_exceeded`, `envelope_from_corsika_propagation_cpu_budget_exceeded`, `candidate_envelope_cpu_budget_exceeded`, `candidate_envelope_threads1_cpu_budget_exceeded`. [metrics.json](performance.air-shower/metrics.json) · [details.json](performance.air-shower/details.json)

![Wall time per stage: reference and candidates on the same shower input.](performance.air-shower/wall-by-stage.png)

*Wall time per stage: reference and candidates on the same shower input.*

Other figures:

![camera-comparison](performance.air-shower/physics/envelope/camera-comparison.png)
![camera-comparison](performance.air-shower/physics/envelope_vs_unthinned/camera-comparison.png)
![camera-comparison](performance.air-shower/physics/unthinned/camera-comparison.png)

<a id="performance-air-shower-proton"></a>

### Air shower (proton): obdeect chain against sim_telarray — FAIL

Case: Air shower, proton input (20 deg, 300-10000 GeV), same chain and budgets as performance.air-shower; blocked until the input exists

Same comparison as the gamma air shower, on a proton shower input. Reference: CORSIKA IACT → sim_telarray with full optics. Candidates: the obdeect chain up to the camera stage, and the envelope pipe (thinned, QE bound). Stage shares are of the candidate's CPU time.

| Quantity | Value |
|---|---|
| Reference wall s | 9.30 |
| Reference CPU s | 8.85 |
| Reference optics CPU s | 1.07 |
| Reference camera CPU s | 7.78 |
| Candidate CPU s, threads=1 | 25.84 |
| Envelope candidate CPU s, threads=1 | 12.23 |
| Unthinned candidate CPU s, threads=1 | 14.09 |
| Trace-stream CPU / reference optics CPU | 9.80 |
| Thinning gain, CPU | 1.16 |
| Piping gain, CPU | 1.86 |
| Envelope / reference wall, threads=1 | 0.90 |
| Envelope / reference CPU, threads=1 | 1.38 |
| Timing unreliable (1 = loaded machine or noisy repetitions) | 1 |

| Stage | Candidate CPU share, threads=1 | Envelope CPU s, threads=1 | Envelope CPU share, threads=1 |
|---|---|---|---|
| from-corsika | 7% | 1.09 | 9% |
| propagation | 14% | 3.33 | 28% |
| trace-stream | 42% | 1.97 | 17% |
| to-simtel | 6% | 0.07 | 1% |
| sim_telarray | 31% | 5.24 | 45% |

| Candidate | s per shower | CPU ns per photon |
|---|---|---|
| reference | 0.21 | 168.72 |
| candidate, threads=1 | 0.60 | 492.84 |
| envelope candidate, threads=1 | 0.19 | 233.18 |

| Gated check | Actual | Limit | Verdict |
|---|---|---|---|
| `reference_completed` | 1 | 1 ± 0 | PASS |
| `candidate_threads1_completed` | 1 | 1 ± 0 | PASS |
| `candidate_threadsdefault_completed` | 1 | 1 ± 0 | PASS |
| `candidate_envelope_completed` | 1 | 1 ± 0 | PASS |
| `candidate_envelope_threads1_completed` | 1 | 1 ± 0 | PASS |
| `candidate_unthinned_threads1_completed` | 1 | 1 ± 0 | PASS |
| `physics_completed` | 1 | 1 ± 0 | PASS |
| `from_corsika_propagation_cpu_budget_exceeded` | 1 | 0 ± 0 | FAIL |
| `envelope_from_corsika_propagation_cpu_budget_exceeded` | 1 | 0 ± 0 | FAIL |
| `candidate_envelope_wall_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `candidate_envelope_cpu_budget_exceeded` | 1 | 0 ± 0 | FAIL |
| `candidate_envelope_threads1_wall_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `candidate_envelope_threads1_cpu_budget_exceeded` | 1 | 0 ± 0 | FAIL |
| `equal_statistics_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `candidate_envelope_camera_only_missing` | 0 | 0 ± 0 | PASS |
| `candidate_envelope_threads1_camera_only_missing` | 0 | 0 ± 0 | PASS |

`performance.air-shower-proton`: 12/16 gated checks pass; failed: `from_corsika_propagation_cpu_budget_exceeded`, `envelope_from_corsika_propagation_cpu_budget_exceeded`, `candidate_envelope_cpu_budget_exceeded`, `candidate_envelope_threads1_cpu_budget_exceeded`. [metrics.json](performance.air-shower-proton/metrics.json) · [details.json](performance.air-shower-proton/details.json)

![Wall time per stage: reference and candidates on the same proton shower input.](performance.air-shower-proton/wall-by-stage.png)

*Wall time per stage: reference and candidates on the same proton shower input.*

Other figures:

![camera-comparison](performance.air-shower-proton/physics/envelope/camera-comparison.png)
![camera-comparison](performance.air-shower-proton/physics/envelope_vs_unthinned/camera-comparison.png)
![camera-comparison](performance.air-shower-proton/physics/unthinned/camera-comparison.png)

<a id="performance-bottleneck"></a>

### Bottleneck of the envelope pipe — PASS

Case: Bottleneck of the envelope pipe (LSTN-01 7.0.0, 20 deg gamma): prefix pipes from-corsika to sim_telarray, marginal wall per member (report only)

Prefix pipes: the pipe is run with its first 1, 2, … 5 members, the last one writing to /dev/null. The wall time a member adds (marginal) is what it costs inside the pipe; the member with the largest one is the bottleneck. CPU s is the member's own CPU time in the full pipe.

| Quantity | Value |
|---|---|
| Bottleneck member | propagation |
| Full pipe wall s (median) | 3.41 |
| Full pipe wall spread over repetitions | 2.0% |

| Member | Marginal wall s | Wall share | CPU s |
|---|---|---|---|
| from-corsika | 0.34 | 10% | 0.37 |
| propagation | 1.26 | 37% | 1.24 |
| trace-stream | 0.54 | 16% | 1.05 |
| to-simtel | 0.07 | 2% | 0.03 |
| sim_telarray | 1.19 | 35% | 2.46 |

| Gated check | Actual | Limit | Verdict |
|---|---|---|---|
| `prefix1_completed` | 1 | 1 ± 0 | PASS |
| `prefix2_completed` | 1 | 1 ± 0 | PASS |
| `prefix3_completed` | 1 | 1 ± 0 | PASS |
| `prefix4_completed` | 1 | 1 ± 0 | PASS |
| `prefix5_completed` | 1 | 1 ± 0 | PASS |

`performance.bottleneck`: 5/5 gated checks pass. [metrics.json](performance.bottleneck/metrics.json) · [details.json](performance.bottleneck/details.json)

![Wall time per member of the envelope pipe, and the bottleneck.](performance.bottleneck/bottleneck.png)

*Wall time per member of the envelope pipe, and the bottleneck.*

<a id="performance-stages"></a>

### Stage timing (wall and CPU) — PASS

Case: Per-stage throughput on the fixed CORSIKA IACT input, air-shower setup (LSTN-01 7.0.0; report only, no performance gates)

Each chain stage on one CORSIKA input, repeated; the table gives the medians, the wall-time spread over the repetitions and each stage's share of the CPU time.

| Quantity | Value |
|---|---|
| Timing unreliable (1 = loaded machine or noisy repetitions) | 1 |
| Host load average before the run | 6.77 |

| Stage | Wall s (median) | CPU s (median) | Photons/s | CPU share | CPU ns per photon | Wall spread |
|---|---|---|---|---|---|---|
| from-corsika | 1.46 | 1.12 | 3.72e+07 | 17% | 21 | 83.0% |
| propagation | 1.65 | 1.27 | 3.29e+07 | 19% | 23 | 24.3% |
| propagation (envelope) | 1.36 | 1.24 | 3.98e+07 | 18% | 23 | 10.1% |
| compile model | 0.52 | 0.41 | – | – | – | 47.7% |
| trace-stream | 3.83 | 3.92 | 9.9e+06 | 58% | 103 | 2.9% |
| to-simtel | 0.54 | 0.50 | 1.01e+08 | 7% | 9 | 54.6% |
| sim_telarray (direct) | 4.29 | 3.98 | 1.26e+07 | 59% | 73 | 20.1% |

| Gated check | Actual | Limit | Verdict |
|---|---|---|---|
| `stages_failed` | 0 | 0 ± 0 | PASS |
| `stages_completed` | 6 | 6 ± 0 | PASS |
| `from_corsika_completed` | 1 | 1 ± 0 | PASS |
| `to_simtel_completed` | 1 | 1 ± 0 | PASS |
| `propagation_completed` | 1 | 1 ± 0 | PASS |
| `propagation_envelope_completed` | 1 | 1 ± 0 | PASS |
| `compile_model_completed` | 1 | 1 ± 0 | PASS |
| `trace_stream_completed` | 1 | 1 ± 0 | PASS |
| `simtel_direct_completed` | 1 | 1 ± 0 | PASS |
| `from_corsika_reps` | 3 | 3 ± 0 | PASS |
| `to_simtel_reps` | 3 | 3 ± 0 | PASS |
| `propagation_reps` | 3 | 3 ± 0 | PASS |
| `propagation_envelope_reps` | 3 | 3 ± 0 | PASS |
| `compile_model_reps` | 3 | 3 ± 0 | PASS |
| `trace_stream_reps` | 3 | 3 ± 0 | PASS |
| `simtel_direct_reps` | 3 | 3 ± 0 | PASS |
| `iact_input_vs_stream_weight_relative_delta` | 1.847e-11 | 0 ± 1e-06 | PASS |
| `from_corsika_to_simtel_weight_relative_delta` | 1.847e-11 | 0 ± 1e-06 | PASS |

`performance.stages`: 18/18 gated checks pass. [metrics.json](performance.stages/metrics.json) · [details.json](performance.stages/details.json)

![Photons per second for each stage.](performance.stages/photons-per-second.png)

*Photons per second for each stage.*

<a id="performance-scaling"></a>

### Scaling: fixed and per-shower cost, trace-stream threads — PASS

Case: Air shower, one LST (LSTN-01 7.0.0, 20 deg gamma): fixed vs per-shower cost of each chain stage (0, N, 2N showers) and trace-stream thread scaling (report only)

Wall time of each stage split into a fixed part and a part per shower. Trace-stream speedup and efficiency for 1, 2, 4 and 8 threads.

| Quantity | Value |
|---|---|
| Cores on the host | 6 |

| Stage | Fixed wall s | Per-shower wall s | Fixed CPU s | Per-shower CPU s | Fixed share at N |
|---|---|---|---|---|---|
| from-corsika | -0.04 | 0.0256 | -0.01 | 0.0187 | -16% |
| propagation | -0.02 | 0.0637 | 0.01 | 0.0528 | -3% |
| compile model | 0.34 | 0 | 0.47 | 0 | 100% |
| trace-stream | 0.24 | 0.156 | 0.55 | 0.14 | 13% |
| to-simtel | -0.00 | 0.0128 | 0.00 | 0.00999 | -3% |

| Trace-stream threads | Wall s | CPU s | Speedup | Efficiency |
|---|---|---|---|---|
| 1 | 1.70 | 1.83 | 1.00 | 100% |
| 2 | 1.13 | 1.75 | 1.50 | 75% |
| 4 | 0.89 | 1.88 | 1.91 | 48% |

| Gated check | Actual | Limit | Verdict |
|---|---|---|---|
| `from_corsika_scaling_completed` | 1 | 1 ± 0 | PASS |
| `propagation_scaling_completed` | 1 | 1 ± 0 | PASS |
| `compile_model_scaling_completed` | 1 | 1 ± 0 | PASS |
| `trace_stream_scaling_completed` | 1 | 1 ± 0 | PASS |
| `to_simtel_scaling_completed` | 1 | 1 ± 0 | PASS |
| `thread_scaling_completed` | 1 | 1 ± 0 | PASS |

`performance.scaling`: 6/6 gated checks pass. [metrics.json](performance.scaling/metrics.json) · [details.json](performance.scaling/details.json)

![Wall time against number of showers per stage, and trace-stream speedup against threads.](performance.scaling/scaling.png)

*Wall time against number of showers per stage, and trace-stream speedup against threads.*

<a id="performance-air-shower-array"></a>

### Array of six telescopes: chain against one sim_telarray run — PASS

Case: Air shower, North array (LSTN-01..04, MSTN-03, MSTN-15; 7.0.0): CORSIKA IACT -> sim_telarray with all six telescopes versus an obdeect tee candidate (one QE-envelope propagation per telescope group, FIFO fan-out to one sim_telarray camera stage per telescope, Bypass_Optics=3); chain time without compiles; per-telescope photoelectron bias; timing budgets

Reference: one sim_telarray run with all telescopes, full optics. Candidate: the obdeect chain per telescope, one pipe per telescope, and the variant with a tee (fan-out).

| Quantity | Value |
|---|---|
| Candidate with tee, wall s | 10.79 |
| Candidate with tee, CPU s | 24.44 |
| Reference wall s | 45.33 |
| Reference CPU s | 43.55 |
| Tee / per-telescope CPU | 0.76 |
| Timing unreliable (1 = loaded machine or noisy repetitions) | 1 |

| Gated check | Actual | Limit | Verdict |
|---|---|---|---|
| `reference_completed` | 1 | 1 ± 0 | PASS |
| `candidate_tee_completed` | 1 | 1 ± 0 | PASS |
| `physics_completed` | 1 | 1 ± 0 | PASS |
| `array_chain_wall_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `array_chain_cpu_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `equal_statistics_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `photon_gate_budget_exceeded` | 0 | 0 ± 0 | PASS |
| `tel1_camera_only_missing` | 0 | 0 ± 0 | PASS |
| `tel2_camera_only_missing` | 0 | 0 ± 0 | PASS |
| `tel3_camera_only_missing` | 0 | 0 ± 0 | PASS |
| `tel4_camera_only_missing` | 0 | 0 ± 0 | PASS |
| `tel5_camera_only_missing` | 0 | 0 ± 0 | PASS |
| `tel6_camera_only_missing` | 0 | 0 ± 0 | PASS |

`performance.air-shower-array`: 13/13 gated checks pass. [metrics.json](performance.air-shower-array/metrics.json) · [details.json](performance.air-shower-array/details.json)

![Wall time per stage, reference and candidate.](performance.air-shower-array/wall-by-stage.png)

*Wall time per stage, reference and candidate.*

Other figures:

![camera-comparison](performance.air-shower-array/physics/candidate_per_telescope_tel1-LSTN-01/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/candidate_per_telescope_tel2-LSTN-02/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/candidate_per_telescope_tel3-LSTN-03/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/candidate_per_telescope_tel4-LSTN-04/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/candidate_per_telescope_tel5-MSTN-03/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/candidate_per_telescope_tel6-MSTN-15/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/tel1-LSTN-01/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/tel2-LSTN-02/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/tel3-LSTN-03/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/tel4-LSTN-04/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/tel5-MSTN-03/camera-comparison.png)
![camera-comparison](performance.air-shower-array/physics/tel6-MSTN-15/camera-comparison.png)

<a id="performance-stream-pipe"></a>

### Stream pipes: throughput — FAIL

Case: Stream pipe throughput, budget >= 1 GB/s: cat | tee and cat | obdeect stream to-simtel on the CORSIKA native stream

The native photon stream of the CORSIKA input read by cat and piped into tee (copy to disk) or into obdeect stream to-simtel. Throughput in decimal GB/s; the budget is 1 GB/s.

| Quantity | Value |
|---|---|
| cat → tee, GB/s | 0.81 |
| cat → obdeect stream to-simtel, GB/s | 0.48 |

| Gated check | Actual | Limit | Verdict |
|---|---|---|---|
| `cat_tee_pipe_completed` | 1 | 1 ± 0 | PASS |
| `to_simtel_pipe_completed` | 1 | 1 ± 0 | PASS |
| `cat_tee_pipe_throughput_budget_exceeded` | 1 | 0 ± 0 | FAIL |
| `to_simtel_pipe_throughput_budget_exceeded` | 1 | 0 ± 0 | FAIL |

`performance.stream-pipe`: 2/4 gated checks pass; failed: `cat_tee_pipe_throughput_budget_exceeded`, `to_simtel_pipe_throughput_budget_exceeded`. [metrics.json](performance.stream-pipe/metrics.json) · [details.json](performance.stream-pipe/details.json)

<a id="performance-flasher-dish"></a>

### Flasher on the dish: chain against LightEmission — PASS

Case: Flasher-dish, wall time budget: LightEmission -> sim_telarray versus obdeect sources -> propagation -> stream -> sim_telarray, matched photon counts (candidate wall <= reference wall)

Reference: LightEmission (ff-1m) → sim_telarray. Candidate: obdeect light source → chain → sim_telarray, one thread.

| Quantity | Value |
|---|---|
| Reference wall s | 0.32 |
| Candidate wall s, threads=1 | 0.22 |
| Candidate / reference wall, threads=1 | 0.70 |

| Stage | Wall s, threads=1 |
|---|---|
| light source | 0.01 |
| propagation | 0.00 |
| to-simtel | 0.01 |
| sim_telarray | 0.20 |

| Gated check | Actual | Limit | Verdict |
|---|---|---|---|
| `reference_completed` | 1 | 1 ± 0 | PASS |
| `candidate_threads1_completed` | 1 | 1 ± 0 | PASS |
| `candidate_threads1_wall_budget_exceeded` | 0 | 0 ± 0 | PASS |

`performance.flasher-dish`: 3/3 gated checks pass. [metrics.json](performance.flasher-dish/metrics.json) · [details.json](performance.flasher-dish/details.json)

![Wall time per stage, reference and candidate.](performance.flasher-dish/wall-by-stage.png)

*Wall time per stage, reference and candidate.*

<a id="performance-point-source"></a>

### Point source: chain against LightEmission — PASS

Case: Point-source, wall time budget: LightEmission -> sim_telarray versus obdeect sources -> propagation -> stream -> sim_telarray, matched photon counts (candidate wall <= reference wall)

Reference: LightEmission → sim_telarray. Candidate: obdeect light source → chain → sim_telarray, one thread.

| Quantity | Value |
|---|---|
| Reference wall s | 0.19 |
| Candidate wall s, threads=1 | 0.15 |
| Candidate / reference wall, threads=1 | 0.78 |

| Stage | Wall s, threads=1 |
|---|---|
| light source | 0.00 |
| propagation | 0.00 |
| to-simtel | 0.00 |
| sim_telarray | 0.15 |

| Gated check | Actual | Limit | Verdict |
|---|---|---|---|
| `reference_completed` | 1 | 1 ± 0 | PASS |
| `candidate_threads1_completed` | 1 | 1 ± 0 | PASS |
| `candidate_threads1_wall_budget_exceeded` | 0 | 0 ± 0 | PASS |

`performance.point-source`: 3/3 gated checks pass. [metrics.json](performance.point-source/metrics.json) · [details.json](performance.point-source/details.json)

![Wall time per stage, reference and candidate.](performance.point-source/wall-by-stage.png)

*Wall time per stage, reference and candidate.*
