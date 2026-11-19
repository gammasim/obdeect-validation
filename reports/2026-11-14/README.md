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

[6. Performance: obdeect chain against sim_telarray](#step-performance) is at the bottom: 3 FAIL, 6 PASS.

| Step | Kind | Case | Headline | Verdict |
|---|---|---|---|---|
| 1. Photon stream format | Format | [LED stream → sim_telarray EventIO, read back with pyeventio](#stream-led-eventio) | photons 60 / 60, block order as CORSIKA: 1 | PASS |
| 1. Photon stream format | Format | [CORSIKA → stream → sim_telarray](#stream-corsika-roundtrip) | p.e. +0.07% ± 0.09%, 0 trigger mismatches | PASS |
| 2. Atmosphere (propagation) | Calculation | [Extinction vs analytic Beer–Lambert](#propagation-extinction-analytic) | largest transmission error 0.0e+00 | PASS |
| 2. Atmosphere (propagation) | Calculation | [Extinction vs sim_telarray tables](#propagation-extinction-simtel-table) | largest error 4.5e-08 | PASS |
| 2. Atmosphere (propagation) | Calculation | [Refraction and timing vs CORSIKA atmo](#propagation-atmosphere-corsika) | displacement 4.4e-04 m, time 1.6e-03 ns | PASS |
| 3. Light sources | Monte Carlo | [LED vs pixled](#sources-led-reference) | photons 60 / 60 | PASS |
| 3. Light sources | Monte Carlo | [Flasher on the dish vs ff-1m](#sources-flasher-dish) | photons 1.0015 ± 0.0008, camera p.e. 1.0013 ± 0.0019 | PASS |
| 3. Light sources | Monte Carlo | [Flasher via secondary mirror vs ff-gct](#sources-flasher-secondary) | photons 1.0002 ± 0.0022, time +0.008 ns, camera p.e. 1.001 ± 0.006 | PASS |
| 3. Light sources | Monte Carlo | [Point source vs xyzls](#sources-point-source) | accepted photons 0.995 ± 0.012, time +0.004 ns, camera p.e. 1.010 ± 0.010 | PASS |
| 3. Light sources | Monte Carlo | [Beam scan vs xyls](#sources-beam-scan) | photons 1.000 ± 0.006 | PASS |
| 3. Light sources | Monte Carlo | [Laser beam geometry vs ls-beam](#sources-laser-beam) | path 131.60 vs 131.65 m | PASS |
| 3. Light sources | Monte Carlo | [Night sky vs analytic and nsbls](#sources-night-sky) | rate error 9.2e-08, area-uniform KS 0.0006, camera p.e. vs nsbls 0.80 (recorded only) | PASS |
| 3. Light sources | Monte Carlo | [Muon track vs fake-muon and analytic](#sources-muon-track) | yield 1.0003 ± 0.0029, time vs analytic -1.1e-05 ns, camera p.e. 0.997 ± 0.009 | PASS |
| 4. Telescope optics | Ray tracing | [Optics lst-7.0.0](#optics-lst-7-0-0) | D80 ≤ 0.53%, area ≤ 0.09%, focal length ≤ 0.018% | PASS |
| 4. Telescope optics | Ray tracing | [Optics mst-7.0.0](#optics-mst-7-0-0) | D80 ≤ 1.32%, area ≤ 0.18%, focal length ≤ 0.059% | PASS |
| 4. Telescope optics | Ray tracing | [Optics sst-7.0.0](#optics-sst-7-0-0) | D80 ≤ 0.13%, area ≤ 0.33%, focal length ≤ 0.010% | PASS |
| 5. Full chain | Monte Carlo (shower chain) | [Near-threshold showers (20–100 GeV), stream round trip](#chain-stream-roundtrip-threshold) | 291 telescope events, p.e. +0.00%, 0 trigger mismatches | PASS |
| 5. Full chain | Monte Carlo (shower chain) | [obdeect atmosphere + sim_telarray optics and camera](#chain-atmosphere-substitution) | p.e. +0.12% ± 0.15%, image χ² 1.12 | PASS |
| 5. Full chain | Monte Carlo (shower chain) | [obdeect atmosphere and optics + sim_telarray camera](#chain-optics-bypass1) | p.e. +0.06% ± 0.19%, image χ² 0.99, 0 trigger mismatches | PASS |

## How to read this report

**What is compared.** Each case compares obdeect with a reference: an analytic expectation, a sim_telarray, LightEmission or CORSIKA program, or both. Signed values are obdeect / reference − 1, ratios are obdeect / reference, ± is one standard error.

**What PASS means.** All gated checks of the case are inside limits fixed before the run (`details.json`). Where a reference program is known to be wrong, the case is gated against the analytic expectation, and the comparison with that program is *recorded only*: shown, but not part of the verdict.

**Kind of test.** Format: the data arrive unchanged (no physics). Calculation: a deterministic formula, compared point by point. Monte Carlo: random photons or showers are generated and compared as distributions, with statistical errors. Ray tracing: photons are traced through the telescope optics.

**Same frame before comparing.** obdeect light sources write photons where they are emitted, the LightEmission programs where they arrive. Before comparing, obdeect photons are moved to the reference plane (flight time added). The time panels show the arrival times as compared.

**Known features of the reference programs.** Some differences in the plots come from the reference, not from obdeect; each affected case has an *Expected differences* box under its figure.

- ff-1m, ff-gct and fake-muon write times in steps of 1/32 ns: their time histograms show spikes or a comb, obdeect times are continuous. Means and widths are compared. Cause (from the sim_telarray source): the programs write bunches in the CORSIKA IACT long format, where the time is a 32-bit float relative to a time offset (about 4e5 ns for a 20° 'shower' from the top of the atmosphere, `iact.c`, `toffset`); a float between 2¹⁸ and 2¹⁹ ns has a step of 2⁻⁵ ns = 1/32 ns. Not a physics effect, and 31 ps is small against the pulses.
- nsbls samples the photon entry point wrongly (radius `1 − sqrt(u)` instead of `sqrt(u)`, rotated azimuth): its entry-point, time and camera distributions differ from the analytic expectation, which obdeect follows.
- fake-muon times deviate slightly from the analytic distribution; obdeect follows it.

**Image χ²/ndf (chain cases).** Camera images of the same event and telescope from the two runs, pixel by pixel: mean over pixels with signal of (a − b)² / (a + b), with a, b the photoelectrons of the pixel in each run. Differences from Poisson fluctuations alone give about 1; systematic shape differences give more (limit 2). Values below 1 mean the two runs share part of their fluctuations (e.g. identical photons in the stream round trip). For light sources the summed images are compared the same way (χ² per pixel of the normalised images).

**KS.** Largest difference of two cumulative distributions, compared with the critical value for the sample sizes (in brackets where shown).

**Optics.** The star is moved off axis along x (E/W) and along y (N/S). The optics are not rotationally symmetric, so off axis the two directions give different values; each figure shows both, and both programs must agree for each direction.

## 1. Photon stream format

**Kind of test:** Format.

The obdeect photon stream must carry CORSIKA photons and light-source photons to sim_telarray without loss: CORSIKA → stream → sim_telarray EventIO, compared with sim_telarray on the original CORSIKA file.

<!-- notes:stream -->

<!-- /notes -->

<a id="stream-led-eventio"></a>

### LED stream → sim_telarray EventIO, read back with pyeventio — PASS

Case: LED light-source stream to sim_telarray EventIO, read by the independent eventio reader

An LED light-source stream is converted to sim_telarray EventIO and read back with pyeventio (Python package `eventio`), a reader written independently of both obdeect and sim_telarray, so the file format is not checked by the code that wrote it.

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

> **Expected differences.** The photons are not bit-identical: the stream stores positions in m as float32, so x, y and emission height come back with a relative error up to 0.08 × 1e-6 (the codec bound); directions, times, wavelengths and bunch sizes are exact. A change this small can still move a photon across a mirror or camera edge in sim_telarray and shift its random-number sequence for the rest of the event, so a few pixels differ (image χ²/ndf 0.79, below 1: the runs share most fluctuations). The mean bias (+0.07% ± 0.09%) is zero within errors and no trigger changes. The mechanism is the expected one, not isolated.

## 2. Atmosphere (propagation)

**Kind of test:** Calculation.

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

Extinction from a sim_telarray optical-depth table against sim_telarray's own interpolation, for wavelengths, heights, zenith angles, telescope elevations and in-telescope path lengths. Two obdeect output modes are checked. Weighted: no photon is dropped, each photon's weight is multiplied by its transmission p, so the weight is compared with sim_telarray's interpolated value point by point (largest error above). Thinned: each photon survives with probability p (weight 1); the surviving count is compared with the expectation, hence the residual in σ.

| Quantity | Value |
|---|---|
| Grid (wavelengths × heights × zeniths) | 9 × 5 × 4 |
| Largest error (weighted mode) | 4.5e-08 |
| Largest error (in-telescope depth) | 4.3e-19 |
| Thinned count residual | 2.77 σ |

`propagation.extinction-simtel-table`: 10/10 gated checks pass. [metrics.json](propagation.extinction-simtel-table/metrics.json)

![Weighted mode (photon weight × transmission, no photons dropped): obdeect transmission against the sim_telarray table.](propagation.extinction-simtel-table/extinction-weighted.png)

*Weighted mode (photon weight × transmission, no photons dropped): obdeect transmission against the sim_telarray table.*

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

**Kind of test:** Monte Carlo.

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

**Kind of test:** Ray tracing.

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

![obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = limit.](optics.lst-7.0.0/optics_relative_difference.png)

*obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = limit.*

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

![obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = limit.](optics.mst-7.0.0/optics_relative_difference.png)

*obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = limit.*

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

![obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = limit.](optics.sst-7.0.0/optics_relative_difference.png)

*obdeect / sim_telarray − 1 against off-axis angle, per offset direction; shaded band = limit.*

Plots written by simtools (one program each): [obdeect/validate_optics_SSTS-01_eff_area.png](optics.sst-7.0.0/obdeect/validate_optics_SSTS-01_eff_area.png) · [obdeect/validate_optics_SSTS-01_eff_flen.png](optics.sst-7.0.0/obdeect/validate_optics_SSTS-01_eff_flen.png) · [obdeect/validate_optics_SSTS-01_psf_cm.png](optics.sst-7.0.0/obdeect/validate_optics_SSTS-01_psf_cm.png) · [obdeect/validate_optics_SSTS-01_psf_deg.png](optics.sst-7.0.0/obdeect/validate_optics_SSTS-01_psf_deg.png) · [sim_telarray/validate_optics_SSTS-01_eff_area.png](optics.sst-7.0.0/sim_telarray/validate_optics_SSTS-01_eff_area.png) · [sim_telarray/validate_optics_SSTS-01_eff_flen.png](optics.sst-7.0.0/sim_telarray/validate_optics_SSTS-01_eff_flen.png) · [sim_telarray/validate_optics_SSTS-01_psf_cm.png](optics.sst-7.0.0/sim_telarray/validate_optics_SSTS-01_psf_cm.png) · [sim_telarray/validate_optics_SSTS-01_psf_deg.png](optics.sst-7.0.0/sim_telarray/validate_optics_SSTS-01_psf_deg.png)

## 5. Full chain

**Kind of test:** Monte Carlo (shower chain).

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

## 6. Performance: obdeect chain against sim_telarray

How fast is the obdeect chain on the same showers as sim_telarray, and where does the time go? Report only, apart from a few time budgets (checks below each case).

**What is compared**

Two ways from the same CORSIKA file to camera images:

```
sim_telarray:   CORSIKA file ─▶ sim_telarray [atmosphere + optics + camera]

obdeect chain:  CORSIKA file ─▶ from-corsika ─▶ propagation ─▶ trace-stream ─▶ to-simtel ─▶ sim_telarray [camera only]
                                (read photons)   (atmosphere)   (optics)        (hand over)   (Bypass_Optics=3)
```

The five obdeect programs are connected by a pipe and run at the same time. The camera is the same sim_telarray code on both sides, so the differences come from how the photons reach the camera.

- **QE pre-selection** ("envelope" in the code, "efficiency bound" in obdeect): most Cherenkov photons are never detected, because the quantum efficiency (QE) of the photodetectors is at most about 40 % and much lower at most wavelengths. sim_telarray traces every photon through the optics and drops most of them only in the camera. With QE pre-selection, obdeect's propagation keeps each photon already at the start, with a probability just above the best QE at its wavelength (QE / 0.73); the camera then applies the remaining factor (0.73). The expected signal is unchanged, which each run checks by comparing the camera images (photoelectrons, image shape, arrival times), and about nine of ten photons are dropped before the optics. The saving in computing time is much smaller than that: without pre-selection the photons stay in CORSIKA's bunches of several photons, and trace-stream itself drops most of them just before tracing (by the mirror reflectivity), so the work per shower changes less than the photon count.
- Same showers, same telescope models (CTAO North, production 7.0.0, 20° zenith), same machine. Every obdeect program runs with one thread unless a table says otherwise.

**How to read the times**

- **CPU time** is the processor time summed over all programs: what a batch farm pays. **Wall time** is the elapsed time. In a pipe the five obdeect programs run at the same time on different cores, so the chain's wall time is shorter than its CPU time; on a single core it would be at least the CPU time.
- Ratios are obdeect chain / sim_telarray: below 1 means obdeect is faster.
- Building the optical model (`ray-tracing: compile`, about 0.3 s per run) is left out of the chain times unless stated; it is listed on its own.
- **Timing unreliable** = 1: the machine was busy (load above 1) or the repetitions differ by more than 5 %. Repetitions alternate between sim_telarray and obdeect, so both see the same conditions; the "fastest repetition" ratios use the fastest run of each, closest to an idle machine.

Host: Apple M3 Pro (6 performance + 6 efficiency cores), podman VM 6 CPUs

<!-- notes:performance -->
**Run (macOS).** 18–19 November 2026: the nine cases, **everything on one core** (`CPUS=0`), so wall and CPU time are comparable and the five programs of a pipe share that core. obdeect branch `native-trace-stream` plus uncommitted changes (`--telescopes` for propagation in the harness, `StdinBuffer` for `propagation run` and `trace-stream`); sim_telarray with `Bypass_Optics=3`. Machine: Apple M3 Pro laptop, podman VM with 6 CPUs, shared with other users (load 1–2.7 before the runs). Every case is flagged "timing unreliable" (load above 1 or repetitions differing by more than 5 %), but the gamma case was measured in all ten suite calls of the day: sim_telarray 14.3–14.7 s CPU, chain 14.1–14.7 s, ratio 0.97–1.01. Final numbers need an idle node.

**In short**

1. On one LST, with the same showers and the same camera images (photoelectron bias −0.04 % ± 0.15 %), the obdeect chain needs **0.99× (gamma) and 0.96× (proton) the CPU time** of sim_telarray and 1.02× and 1.00× the wall time. The CPU budget of one thread per program passes, the wall budget (≤ 1.00) misses by 2 % for gammas. It is a tie, not a win: single core, QE pre-selection, same camera.
2. **Three quarters of the CPU is the sim_telarray camera simulation** (10.6 of 14.1 s for gammas), the same code on both sides. It receives fewer photons in the obdeect chain than in the sim_telarray run (QE pre-selection): 2.2 s less. That pays for the rest of the chain: from-corsika, propagation, trace-stream and to-simtel take 3.6 s, sim_telarray's optics with the atmosphere 1.4 s.
3. **Reading and atmosphere (from-corsika + propagation) are 20 % (gamma) and 30 % (proton) of the sim_telarray CPU**; the budget is 10 %. This is the part to optimise next.
4. **The optics are not the problem.** trace-stream traces the preselected photons in 0.83 s (sim_telarray's optics, with the atmosphere: 1.4 s).
5. **Six telescopes: the chain costs the same CPU per telescope as for one LST**, 13 s on macOS (33 s on Linux). The ratio to sim_telarray (0.43× macOS, 0.59× Linux) is lower only because sim_telarray costs more per telescope in the six-telescope run (30 s against 14 s on macOS, 56 s against 34 s on Linux). This is not an obdeect gain, and the cause is not understood (see the case). Do not quote it as the expected gain for an array.
6. **Checks.** Failing: reading + atmosphere ≤ 10 %, the wall budget and the all-cores CPU budget for gammas, the wall budget and the all-cores budgets for protons, and the 1 GB/s pipe budget, which a plain copy does not reach here either. All physics checks pass.


#### Linux and macOS compared

The same nine cases, run on the Linux machine (apptainer, 4 cases at once, one core each) and on the macOS laptop (podman, one case at a time, one core). Both use the same harness, the same sim_telarray commit (mode 3) and the same obdeect checkout. The CORSIKA showers differ (see **Differences** below), so compare ratios and ns per photon, not seconds.

| | macOS (2026-11-18) | Linux (2026-11-19) |
|---|---|---|
| CPU | Apple M3 Pro (6 P + 6 E cores) | Intel Xeon Gold 6130, 2.1 GHz (64 threads, shared) |
| Cores used per case | 1 | 1 (up to 4 cases at the same time, each pinned to its own core) |
| Container | podman VM, 6 CPUs | apptainer, `simtools-prod_20261118-v78050` |
| CORSIKA | 7.8010, QGSJET-II | 7.8050, QGSJET-III |
| Load before the air-shower case | 1.36 | 2.78 |

**In short**

1. **Same result on both machines.** The obdeect chain needs 1.00× (gamma) and 0.96× (proton) the CPU of sim_telarray on Linux, against 0.99× and 0.96× on macOS. Reading + atmosphere are 18 % and 31 % of the sim_telarray CPU (macOS 20 % and 30 %). The wall budget (≤ 1.00) still fails on Linux, by 10 % for gammas and 6 % for protons (macOS 3 % and 2 %).
2. **The Linux node is about 2× slower per core** than the M3 Pro (gamma 431 against 211 ns per photon, proton 1.6×), for every program alike. It does not change which program dominates: the sim_telarray camera stays at about 75 % of the CPU.
3. **Six telescopes: the chain costs the same per telescope as for one LST** (33 s on Linux, 34 s for one LST). The ratio to sim_telarray, 0.59× on Linux and 0.43× on macOS, comes from sim_telarray costing 1.6× (Linux) and 2.1× (macOS) more per telescope in the six-telescope run. Not understood; it is on the sim_telarray side, so the ratio is not the expected gain for an array.
4. **Pipes are slower in wall time on Linux** (0.26 GB/s against 0.68 GB/s for `cat | tee`) at the same CPU time: the 1.5 GB input comes from the disk of a shared node, not from a warm page cache. Both machines miss the 1 GB/s budget.

![obdeect chain and the single programs, macOS against Linux.](linux-vs-macos.png)

**Is the obdeect chain faster than sim_telarray?** Ratio chain / sim_telarray, below 1 = obdeect faster. The two machines agree.

| Case | Ratio | macOS | Linux |
|---|---|---|---|
| Gamma, one LST | CPU | 0.99 | 1.00 |
| | wall | 1.02 | 1.10 |
| | reading + atmosphere, share of sim_telarray CPU (budget 10 %) | 20 % | 18 % |
| Proton, one LST | CPU | 0.96 | 0.96 |
| | wall | 1.00 | 1.05 |
| | reading + atmosphere, share of sim_telarray CPU (budget 10 %) | 30 % | 31 % |
| Six telescopes | CPU | 0.43 | 0.59 |
| | wall | 0.47 | 0.75 |
| Flasher on the dish | wall | 0.83 | 0.81 |
| Distant point source | wall | 0.74 | 0.76 |

**How fast is each machine?** CPU time per photon that reaches the telescope (gamma case, one LST, one core). Lower is faster.

| | macOS | Linux | Linux / macOS |
|---|---|---|---|
| Gamma: sim_telarray, ns per photon | 211 | 431 | 2.0× |
| Gamma: obdeect chain, ns per photon | 210 | 432 | 2.1× |
| Proton: sim_telarray, ns per photon | 151 | 240 | 1.6× |
| Proton: obdeect chain, ns per photon | 148 | 231 | 1.6× |

**Each program on its own** (the 100-shower gamma file, CPU ns per photon of its input; lower is faster):

| Program | macOS | Linux | Linux / macOS |
|---|---|---|---|
| from-corsika | 5.6 | 11.0 | 2.0× |
| propagation (QE pre-selection) | 6.9 | 13.7 | 2.0× |
| trace-stream | 144.2 | 305.2 | 2.1× |
| to-simtel | 3.6 | 7.9 | 2.2× |
| sim_telarray, full run | 73.4 | 143.4 | 2.0× |

**Photon stream through a pipe** (1 GB/s budget; CPU time is similar, wall time is not):

| | macOS | Linux |
|---|---|---|
| `cat \| tee`: GB/s, wall s, CPU s | 0.68, 1.9, 1.7 | 0.26, 5.8, 1.8 |
| `cat \| to-simtel`: GB/s, wall s, CPU s | 0.87, 1.5, 1.0 | 0.28, 5.5, 2.1 |

**Differences between the two runs, read before comparing seconds**

- **Other showers.** The image has no CORSIKA 7.8010 / QGSJET-II. The Linux inputs were made from the same cards with the image's CORSIKA 7.8050 (`corsika_qgs3_urqmd_flat`, QGSJET-III), so the showers differ: 17 % more photons at the telescope for gammas, 2.1× more for protons. Compare ratios and ns per photon. The CORSIKA logs are in `linux/`.
- **Four cases at once.** Each case ran on its own core, four at a time (cores 10–13), on a node shared with other users (load 3–23 during the runs). Every case is flagged "timing unreliable". The ratios agree with macOS to within 1–10 %; a quiet-node run is still needed for final numbers.
- **Build.** sim_telarray (mode 3, commit b8f2494, hessio 6045783) and obdeect (checkout, 0.9.3.dev7, branch `stdinbuffer`) were built with the image's compilers (GCC 11) under apptainer, not in the podman images. simtools is the workspace checkout.
- **Checks.** The same ones fail as on macOS (reading + atmosphere ≤ 10 %, the wall budget, the all-cores budgets, 1 GB/s pipe). All physics checks pass.

Data: [linux/](linux/) (metrics, details and plots of each case). The macOS data are in the `performance.*` folders.
<!-- /notes -->

| Case | Headline | Verdict |
|---|---|---|
| [Gamma showers, one LST: obdeect chain against sim_telarray](#performance-air-shower) | obdeect / sim_telarray: CPU 0.99, wall 1.02, p.e. bias -0.0% | FAIL |
| [Proton showers, one LST: obdeect chain against sim_telarray](#performance-air-shower-proton) | obdeect / sim_telarray: CPU 0.96, wall 1.00, p.e. bias -0.1% | FAIL |
| [Which program limits the obdeect pipe?](#performance-bottleneck) | bottleneck: sim_telarray, pipe wall 15.2 s | PASS |
| [Each program on its own](#performance-stages) | trace-stream 144 ns/photon, sim_telarray 73 ns/photon (CPU) | PASS |
| [Start-up cost and threads](#performance-scaling) | fixed cost: compile 0.19 s | PASS |
| [Array of six telescopes (4 LST, 2 MST)](#performance-air-shower-array) | obdeect / sim_telarray: CPU 0.43, wall 0.47 | PASS |
| [How fast photons pass between programs](#performance-stream-pipe) | to-simtel 0.87 GB/s, plain copy 0.68 GB/s | FAIL |
| [Flasher on the dish (LST calibration light)](#performance-flasher-dish) | obdeect / sim_telarray side, wall 0.83 | PASS |
| [Distant point source of light](#performance-point-source) | obdeect / sim_telarray side, wall 0.74 | PASS |

<a id="performance-air-shower"></a>

### Gamma showers, one LST: obdeect chain against sim_telarray — FAIL

100 gamma showers (100 GeV–3 TeV, 20° zenith) on LSTN-01. sim_telarray: CORSIKA file → sim_telarray (atmosphere, optics, camera). obdeect chain: one pipe with QE pre-selection, one thread per program. For comparison also timed: the same pipe without QE pre-selection, and the obdeect programs one after another with files in between (no pipe).

| Quantity | Value |
|---|---|
| sim_telarray: wall s, CPU s | 15.09, 14.28 |
| &nbsp;&nbsp;of which optics, CPU s (estimate: full run minus camera-only run) | 1.44 |
| &nbsp;&nbsp;of which camera, CPU s | 12.83 |
| obdeect chain with QE pre-selection: wall s, CPU s | 15.37, 14.08 |
| obdeect chain without QE pre-selection: wall s, CPU s | 16.85, 15.47 |
| obdeect programs one after another (files, no pipe): wall s, CPU s | 24.85, 23.22 |
| ray-tracing: compile, once per run: wall s, CPU s | 0.20, 0.12 |
| **obdeect chain / sim_telarray, CPU** | **0.99** |
| **obdeect chain / sim_telarray, wall** | **1.02** |
| obdeect chain / sim_telarray, CPU, fastest repetitions (compile included) | 1.01 |
| Saving from QE pre-selection (CPU without / with) | 1.10 × |
| Saving from the pipe (CPU of programs one after another / pipe) | 1.49 × |
| Optics: trace-stream CPU without QE pre-selection / sim_telarray optics CPU | 4.71 |
| Timing unreliable | 1 |

| Where the obdeect CPU time goes | With QE pre-selection, CPU s | Share | Without, CPU s | One after another (files), CPU s |
|---|---|---|---|---|
| stream: from-corsika (read the CORSIKA photons) | 1.17 | 8% | 1.11 | 1.26 |
| propagation: run (atmosphere) | 1.51 | 11% | 1.24 | 1.36 |
| ray-tracing: trace-stream (telescope optics) | 0.83 | 6% | 2.22 | 6.79 |
| stream: to-simtel (hand the photons to sim_telarray) | 0.07 | 1% | 0.17 | 0.86 |
| sim_telarray, camera only | 10.59 | 75% | 10.71 | 12.83 |

| Per shower | Wall s per shower | CPU s per shower | CPU ns per photon at the telescope |
|---|---|---|---|
| sim_telarray | 0.151 | 0.143 | 211 |
| obdeect chain with QE pre-selection | 0.156 | 0.142 | 210 |
| obdeect chain without QE pre-selection | 0.170 | 0.156 | 230 |
| obdeect programs one after another | 0.248 | 0.232 | 343 |

| Check | Measured | Result |
|---|---|---|
| Reading + atmosphere (from-corsika + propagation) CPU ≤ 10 % of sim_telarray CPU, programs one after another | 18% | FAIL |
| Reading + atmosphere CPU ≤ 10 % of sim_telarray CPU, in the pipe with QE pre-selection | 20% | FAIL |
| obdeect chain wall ≤ sim_telarray wall (trace-stream on all cores) | 1.04 × sim_telarray | FAIL |
| obdeect chain CPU ≤ sim_telarray CPU (trace-stream on all cores) | 1.00 × sim_telarray | FAIL |
| obdeect chain wall ≤ sim_telarray wall (one thread per program) | 1.02 × sim_telarray | FAIL |
| obdeect chain CPU ≤ sim_telarray CPU (one thread per program) | 0.99 × sim_telarray | PASS |
| Same physics: photoelectron bias of the obdeect chain within its limit | -0.04% ± 0.15% | PASS |

All runs completed with their repetitions (9 checks).

<!-- notes:performance-air-shower -->
**Why FAIL.** Not the physics (photoelectron bias within its limit) and not by much: the chain needs 0.99× the CPU and 1.02× the wall time of sim_telarray. Two budgets set in advance fail: reading + atmosphere ≤ 10 % of the sim_telarray CPU (measured 20 %), and the chain's wall time ≤ sim_telarray's (1.02×, one thread per program). The variant with trace-stream on all cores fails both budgets by 0 to 4 %, because on one core the extra threads only add overhead. These are efficiency targets, not correctness checks.

**Where the CPU goes.** sim_telarray takes 14.28 s: 1.44 s optics with the atmosphere, 12.83 s camera (estimate). The obdeect chain takes 14.08 s:

| Program | CPU s |
|---|---|
| sim_telarray camera | 10.59 |
| propagation | 1.51 |
| from-corsika | 1.17 |
| trace-stream | 0.83 |
| to-simtel | 0.07 |

The camera costs 2.2 s less than in sim_telarray's own run, because it receives only the preselected photons. This pays for reading, the atmosphere and the optics (3.6 s against 1.4 s in sim_telarray).

**QE pre-selection saves 10 %** (15.47 s without it). It drops about nine of ten photons. Without it, the photons stay in CORSIKA's bunches, and trace-stream drops most of them itself before tracing (2.22 s against 0.83 s with it); the extra draw in propagation costs 0.27 s.

**Files.** The obdeect chain writes no intermediate files: its programs are connected by pipes, and only the last program writes (sim_telarray's output, in a container-local work directory, not on the results mount). The "one after another" column is a baseline only, to show what a file-based run would cost: 1.6× the CPU (23.2 s), with 350 MB of intermediate files in the same local directory and a camera that receives photon bunches that were not preselected. It is not a mode of the chain.

**Showers.** 100 CORSIKA showers (seeds 1017…), each used at 5 core positions: 100 telescope events. sim_telarray and the chain read the same file, so the showers are identical within the case. The input file is not hashed in the report (its SHA-256 is in `fixtures/inputs.json`).

**Repeatability.** The three repetitions agree to 1–2 % in CPU time (wall spread up to 10 % for one variant). Across the ten suite calls of the day the ratio was 0.97–1.01 (sim_telarray 14.3–14.7 s CPU, chain 14.1–14.7 s).
<!-- /notes -->

Data: `performance.air-shower` · [metrics.json](performance.air-shower/metrics.json) · [details.json](performance.air-shower/details.json)

![Wall time per program: sim_telarray and the obdeect variants on the same showers.](performance.air-shower/wall-by-stage.png)

*Wall time per program: sim_telarray and the obdeect variants on the same showers.*

Other figures:

![camera-comparison](performance.air-shower/physics/envelope/camera-comparison.png)
![camera-comparison](performance.air-shower/physics/envelope_vs_unthinned/camera-comparison.png)
![camera-comparison](performance.air-shower/physics/unthinned/camera-comparison.png)

<a id="performance-air-shower-proton"></a>

### Proton showers, one LST: obdeect chain against sim_telarray — FAIL

45 proton showers (300 GeV–10 TeV, 20° zenith) on LSTN-01. sim_telarray: CORSIKA file → sim_telarray (atmosphere, optics, camera). obdeect chain: one pipe with QE pre-selection, one thread per program. For comparison also timed: the same pipe without QE pre-selection, and the obdeect programs one after another with files in between (no pipe).

| Quantity | Value |
|---|---|
| sim_telarray: wall s, CPU s | 8.44, 7.93 |
| &nbsp;&nbsp;of which optics, CPU s (estimate: full run minus camera-only run) | 1.09 |
| &nbsp;&nbsp;of which camera, CPU s | 6.84 |
| obdeect chain with QE pre-selection: wall s, CPU s | 8.41, 7.64 |
| obdeect chain without QE pre-selection: wall s, CPU s | 9.86, 8.98 |
| obdeect programs one after another (files, no pipe): wall s, CPU s | 16.69, 15.58 |
| ray-tracing: compile, once per run: wall s, CPU s | 0.20, 0.11 |
| **obdeect chain / sim_telarray, CPU** | **0.96** |
| **obdeect chain / sim_telarray, wall** | **1.00** |
| obdeect chain / sim_telarray, CPU, fastest repetitions (compile included) | 0.99 |
| Saving from QE pre-selection (CPU without / with) | 1.18 × |
| Saving from the pipe (CPU of programs one after another / pipe) | 1.72 × |
| Optics: trace-stream CPU without QE pre-selection / sim_telarray optics CPU | 5.12 |
| Timing unreliable | 1 |

| Where the obdeect CPU time goes | With QE pre-selection, CPU s | Share | Without, CPU s | One after another (files), CPU s |
|---|---|---|---|---|
| stream: from-corsika (read the CORSIKA photons) | 0.96 | 13% | 0.95 | 1.23 |
| propagation: run (atmosphere) | 1.15 | 15% | 0.96 | 1.05 |
| ray-tracing: trace-stream (telescope optics) | 0.67 | 9% | 1.89 | 5.57 |
| stream: to-simtel (hand the photons to sim_telarray) | 0.06 | 1% | 0.16 | 0.78 |
| sim_telarray, camera only | 4.77 | 62% | 5.02 | 6.84 |

| Per shower | Wall s per shower | CPU s per shower | CPU ns per photon at the telescope |
|---|---|---|---|
| sim_telarray | 0.188 | 0.176 | 151 |
| obdeect chain with QE pre-selection | 0.191 | 0.172 | 148 |
| obdeect chain without QE pre-selection | 0.223 | 0.202 | 173 |
| obdeect programs one after another | 0.371 | 0.346 | 297 |

| Check | Measured | Result |
|---|---|---|
| Reading + atmosphere (from-corsika + propagation) CPU ≤ 10 % of sim_telarray CPU, programs one after another | 29% | FAIL |
| Reading + atmosphere CPU ≤ 10 % of sim_telarray CPU, in the pipe with QE pre-selection | 30% | FAIL |
| obdeect chain wall ≤ sim_telarray wall (trace-stream on all cores) | 1.05 × sim_telarray | FAIL |
| obdeect chain CPU ≤ sim_telarray CPU (trace-stream on all cores) | 1.01 × sim_telarray | FAIL |
| obdeect chain wall ≤ sim_telarray wall (one thread per program) | 1.00 × sim_telarray | PASS |
| obdeect chain CPU ≤ sim_telarray CPU (one thread per program) | 0.96 × sim_telarray | PASS |
| Same physics: photoelectron bias of the obdeect chain within its limit | -0.09% ± 0.14% | PASS |

All runs completed with their repetitions (9 checks).

<!-- notes:performance-air-shower-proton -->
**Same picture as for gammas: 0.96× the CPU, 1.00× the wall time.** sim_telarray takes 7.93 s, of which 1.09 s is optics. The obdeect chain takes 7.64 s:

| Program | CPU s |
|---|---|
| sim_telarray camera | 4.77 |
| propagation | 1.15 |
| from-corsika | 0.96 |
| trace-stream | 0.67 |
| to-simtel | 0.06 |

Reading and the atmosphere are 30 % of sim_telarray's CPU (check), more than for gammas (20 %): the camera is cheaper for protons, reading and propagation are not. The case reads the 45-shower proton file, not a 100-shower one.
<!-- /notes -->

Data: `performance.air-shower-proton` · [metrics.json](performance.air-shower-proton/metrics.json) · [details.json](performance.air-shower-proton/details.json)

![Wall time per program: sim_telarray and the obdeect variants on the same showers.](performance.air-shower-proton/wall-by-stage.png)

*Wall time per program: sim_telarray and the obdeect variants on the same showers.*

Other figures:

![camera-comparison](performance.air-shower-proton/physics/envelope/camera-comparison.png)
![camera-comparison](performance.air-shower-proton/physics/envelope_vs_unthinned/camera-comparison.png)
![camera-comparison](performance.air-shower-proton/physics/unthinned/camera-comparison.png)

<a id="performance-bottleneck"></a>

### Which program limits the obdeect pipe? — PASS

The obdeect pipe with QE pre-selection (100 gamma showers, LSTN-01) is built up one program at a time: from-corsika alone, then from-corsika → propagation, and so on, the last program writing its output nowhere. The wall time a program adds is what it costs inside the pipe; the program that adds most limits the whole chain (the bottleneck). CPU s is the program's own CPU time in the full pipe.

| Quantity | Value |
|---|---|
| Bottleneck | sim_telarray |
| Full pipe, wall s | 15.23 |
| Spread of the full pipe over repetitions | 0.8% |

| Program | Wall s it adds | Share of the pipe wall | Its CPU s |
|---|---|---|---|
| stream: from-corsika (read the CORSIKA photons) | 0.99 | 6% | 1.08 |
| propagation: run (atmosphere) | 1.69 | 11% | 1.50 |
| ray-tracing: trace-stream (telescope optics) | 0.88 | 6% | 0.82 |
| stream: to-simtel (hand the photons to sim_telarray) | 0.15 | 1% | 0.07 |
| sim_telarray, camera only | 11.52 | 76% | 10.50 |

All runs completed with their repetitions (5 checks).

<!-- notes:performance-bottleneck -->
**On one core this table shows cost shares, not a bottleneck**: the five programs take turns on the core, so the wall time a program adds is about its CPU time. The sim_telarray camera takes 11.5 of the 15.2 s (76 %), propagation 1.7 s, from-corsika 1.0 s, trace-stream 0.9 s.

On several cores the programs of the pipe run at the same time, and the slowest one sets the wall time: the sim_telarray camera (10.5 s CPU) against 3.5 s for all other programs together. So faster propagation or reading would not shorten the wall time of a multi-core pipe, but it saves CPU time in full, which is what a batch farm pays. This statement is derived from the single-core numbers; a multi-core run was not made for this report.
<!-- /notes -->

Data: `performance.bottleneck` · [metrics.json](performance.bottleneck/metrics.json) · [details.json](performance.bottleneck/details.json)

![Wall time each program adds to the pipe (dark) and its CPU time (light).](performance.bottleneck/bottleneck.png)

*Wall time each program adds to the pipe (dark) and its CPU time (light).*

<a id="performance-stages"></a>

### Each program on its own — PASS

Each obdeect program run alone on the 100 gamma showers (input from a file, output to a file), repeated, with sim_telarray (optics and camera) on the same showers for scale. trace-stream here traces the photons without QE pre-selection, about ten times more than in the production pipe. CPU ns per photon: per photon at the telescope (before the atmosphere).

| Quantity | Value |
|---|---|
| Timing unreliable | 1 |
| Load average before the run | 1.66 |

| Program | Wall s | CPU s | Photons/s | CPU ns per photon | Spread over repetitions |
|---|---|---|---|---|---|
| stream: from-corsika (read the CORSIKA photons) | 1.68 | 1.12 | 1.19e+08 | 6 | 32% |
| propagation: run (atmosphere) | 1.32 | 1.21 | 1.51e+08 | 6 | 30% |
| propagation: run with QE pre-selection | 1.45 | 1.38 | 1.38e+08 | 7 | 6% |
| ray-tracing: compile (build the optical model) | 0.18 | 0.10 | – | – | 9% |
| ray-tracing: trace-stream (telescope optics) | 6.78 | 6.50 | 6.64e+06 | 144 | 6% |
| stream: to-simtel (hand the photons to sim_telarray) | 1.05 | 0.71 | 1.9e+08 | 4 | 30% |
| sim_telarray, optics and camera (production) | 15.52 | 14.64 | 1.28e+07 | 73 | 2% |

| Check | Measured | Result |
|---|---|---|
| Photon count unchanged from CORSIKA file to obdeect stream | 1.2e-10 | PASS |
| Photon count unchanged through from-corsika → to-simtel | 1.2e-10 | PASS |

All runs completed with their repetitions (16 checks).

<!-- notes:performance-stages -->
CPU time of each program on its own (100 gamma showers, one core):

| Program | CPU s |
|---|---|
| from-corsika | 1.12 |
| propagation | 1.21 (1.38 with QE pre-selection) |
| trace-stream (all photons, no pre-selection) | 6.50 |
| to-simtel | 0.71 |
| sim_telarray, optics and camera | 14.64 |

The ns per photon are not comparable between rows: each program is divided by the photons it receives (trace-stream gets only the photons at LSTN-01, from-corsika and sim_telarray the whole file). trace-stream is the most expensive obdeect program here because this case feeds it photons that were not preselected; in the production pipe it sees only the preselected ones (0.83 s). The short programs (from-corsika, to-simtel, propagation) differ by up to 32 % between repetitions on this machine; the pipe cases, which run them for longer, agree to 1–2 %.
<!-- /notes -->

Data: `performance.stages` · [metrics.json](performance.stages/metrics.json) · [details.json](performance.stages/details.json)

![Photons per second for each program.](performance.stages/photons-per-second.png)

*Photons per second for each program.*

<a id="performance-scaling"></a>

### Start-up cost and threads — PASS

Each obdeect program timed on 0, 10 and 20 gamma showers. A straight line through the three points splits its time into a fixed part (start-up, reading tables, building the model) and a part per shower; small negative fixed parts are noise around zero. Below: ray-tracing: trace-stream with 1, 2, 4 (and 8 if the machine has them) threads on the same photons.

| Quantity | Value |
|---|---|
| Cores | 1 |

| Program | Fixed wall s | Wall s per shower | Fixed CPU s | CPU s per shower | Fixed share for 10 showers |
|---|---|---|---|---|---|
| stream: from-corsika (read the CORSIKA photons) | 0.02 | 0.009 | 0.02 | 0.009 | 15% |
| propagation: run (atmosphere) | 0.01 | 0.014 | 0.01 | 0.013 | 9% |
| ray-tracing: compile (build the optical model) | 0.19 | 0.000 | 0.10 | 0.000 | 100% |
| ray-tracing: trace-stream (telescope optics) | 0.06 | 0.067 | 0.06 | 0.063 | 8% |
| stream: to-simtel (hand the photons to sim_telarray) | 0.03 | 0.005 | 0.01 | 0.005 | 33% |

| trace-stream threads | Wall s | CPU s | Speedup | Efficiency |
|---|---|---|---|---|
| 1 | 0.84 | 0.80 | 1.00 | 100% |

All runs completed with their repetitions (6 checks).

<!-- notes:performance-scaling -->
**Fixed cost per run.** Only one program has a noticeable fixed cost: building the optical model (`ray-tracing: compile`), 0.19 s wall and 0.10 s CPU on every run. It is still Python. Everything else starts within a few milliseconds. With hundreds of showers per production run, the fixed cost does not matter. The compiled model could also be built once and reused.

**Per shower** (CPU): from-corsika 0.009 s, propagation 0.013 s, trace-stream 0.063 s, to-simtel 0.005 s.

**Threads.** Not measured in this run: it was pinned to one core, so only 1 thread was timed (speedup 1.00 by definition). The thread scaling needs a multi-core run.
<!-- /notes -->

Data: `performance.scaling` · [metrics.json](performance.scaling/metrics.json) · [details.json](performance.scaling/details.json)

![Left: wall time against the number of showers, with the line fits. Right: trace-stream speedup against threads.](performance.scaling/scaling.png)

*Left: wall time against the number of showers, with the line fits. Right: trace-stream speedup against threads.*

<a id="performance-air-shower-array"></a>

### Array of six telescopes (4 LST, 2 MST) — PASS

100 gamma showers on LSTN-01 to -04, MSTN-03 and MSTN-15. sim_telarray: one run that simulates all six telescopes (atmosphere, optics, cameras), single-threaded. obdeect chain: telescopes with the same QE share one from-corsika → propagation (QE pre-selection); a `tee` copies that photon stream to one trace-stream → to-simtel → sim_telarray (camera only) branch per telescope. The branches are separate processes that run at the same time, so the CPU ratio is the fair comparison; the wall time profits from the parallel branches. Not compared: the array trigger (each telescope has its own sim_telarray run).

| Quantity | Value |
|---|---|
| sim_telarray, all six telescopes: wall s, CPU s | 189.0, 181.8 |
| obdeect chain with tee: wall s, CPU s | 88.2, 77.8 |
| **obdeect chain / sim_telarray, CPU** | **0.43** |
| **obdeect chain / sim_telarray, wall** | **0.47** |
| obdeect chain / sim_telarray, CPU, fastest repetitions | 0.43 |
| Shared propagation (tee) / one full pipe per telescope, CPU | 0.98 |
| Timing unreliable | 1 |

| Telescope | Photoelectron bias obdeect vs sim_telarray | ± (1 σ) | Images compared | sim_telarray camera CPU s | trace-stream CPU s |
|---|---|---|---|---|---|
| LSTN-01 | -0.05% | 0.21% | 484 | 10.20 | 0.53 |
| LSTN-02 | +0.37% | 0.25% | 486 | 10.66 | 0.56 |
| LSTN-03 | +0.11% | 0.24% | 469 | 10.53 | 0.55 |
| LSTN-04 | +0.66% | 0.24% | 478 | 10.59 | 0.53 |
| MSTN-03 | -0.16% | 0.37% | 276 | 11.29 | 0.23 |
| MSTN-15 | -0.93% | 0.32% | 433 | 12.95 | 0.28 |

| Check | Measured | Result |
|---|---|---|
| obdeect chain wall ≤ sim_telarray wall | 0.47 × sim_telarray | PASS |
| obdeect chain CPU ≤ sim_telarray CPU | 0.43 × sim_telarray | PASS |
| Same physics: photoelectron bias within its limit in all six telescopes | 0 | PASS |
| Photon counts consistent through the chain in all six telescopes | 0 | PASS |

All runs completed with their repetitions (9 checks).

<!-- notes:performance-air-shower-array -->
**The chain needs 0.43× the CPU of sim_telarray (77.8 s against 181.8 s) on macOS, but that is not an obdeect gain.** On Linux the ratio is 0.59× (200 s against 337 s). Where the obdeect time goes (macOS, CPU s): six camera stages 66 (10.2–13.0 each), from-corsika and propagation, shared by all telescopes, about 8, trace-stream 3, tee 0.4.

**Per telescope**, in CPU seconds:

| | One LST | Six telescopes |
|---|---|---|
| sim_telarray | 14.3 | 30.3 (181.8 / 6) |
| obdeect chain | 14.1 | 13.0 (77.8 / 6) |

The obdeect cost per telescope is the same as for one LST. Sharing the reading and the atmosphere between telescopes saves little: the shared chain needs 0.98× the CPU of one full pipe per telescope. What differs is sim_telarray: **2.1× more CPU per telescope** (macOS; 1.6× on Linux, 56 s against 34 s) in the six-telescope run than in the one-LST run. If its cameras cost what the obdeect-side cameras cost (66 s in total on macOS), its optics and atmosphere would take about 115 s, 19 s per telescope against 1.4 s for one LST. **This is not understood.** It is not the number of photons per telescope: the array input puts fewer photons on each telescope (mean 1.0e5 per telescope and event) than the single-LST input (1.8e5). Open candidates: per-telescope overhead of a sim_telarray run with six telescopes (trigger, array configuration, the MST models), or the larger core scatter (250 m against 100 m). A sim_telarray run for one telescope of this array on this input, and runs with 1, 2, 4 and 6 identical LSTs, would separate the causes. Until then, the safe statement is: the obdeect chain costs about the same CPU per telescope as sim_telarray for one LST, and sharing the stream between telescopes gives no extra gain. Do not quote 0.43 or 0.59 as the expected gain for an array.

**The showers of the two cases differ.** Single-LST input: seeds 1017…, core scatter 100 m, three telescopes of 12 m radius. Array input: seeds 1617…, scatter 250 m, six telescopes. Within a case, sim_telarray and the obdeect chain read the same file (100 telescope events), so the showers are identical there by construction. Times of different cases are not strictly like for like.

**Physics.** The photoelectron bias is within its limit for all six telescopes (−0.9 % to +0.7 %).

<!-- /notes -->

Data: `performance.air-shower-array` · [metrics.json](performance.air-shower-array/metrics.json) · [details.json](performance.air-shower-array/details.json)

![Wall time per program, sim_telarray and obdeect.](performance.air-shower-array/wall-by-stage.png)

*Wall time per program, sim_telarray and obdeect.*

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

### How fast photons pass between programs — FAIL

The photon stream of the 100 gamma showers (about 1.3 GB) read with `cat` and written to disk by `tee` (a plain copy: the ceiling of this machine), and converted by `stream: to-simtel` (the last obdeect step before sim_telarray). Throughput in GB/s (10^9 bytes per second); the budget is 1 GB/s.

| Quantity | Value |
|---|---|
| Plain copy (cat → tee), GB/s | 0.68 |
| cat → stream: to-simtel, GB/s | 0.87 |

| Check | Measured | Result |
|---|---|---|
| Plain copy ≥ 1 GB/s | 0.68 GB/s | FAIL |
| stream: to-simtel ≥ 1 GB/s | 0.87 GB/s | FAIL |

All runs completed with their repetitions (2 checks).

<!-- notes:performance-stream-pipe -->
A plain copy reaches only 0.68 GB/s in this VM (the podman VM's file system), so the 1 GB/s budget cannot be met on this machine. to-simtel reaches 0.87 GB/s. Its CPU time, 1.0 s for 1.3 GB (about 1.3 GB/s), shows that the conversion itself would keep up; the rest is waiting for the pipe or the disk. This case needs a Linux node with a faster disk.
<!-- /notes -->

Data: `performance.stream-pipe` · [metrics.json](performance.stream-pipe/metrics.json) · [details.json](performance.stream-pipe/details.json)

<a id="performance-flasher-dish"></a>

### Flasher on the dish (LST calibration light) — PASS

2 flashes of 20 million photons each from a flasher on the dish; about 56 000 photons reach the camera. sim_telarray side: LightEmission `ff-1m` writes an uncompressed photon file (as the obdeect stream), then sim_telarray simulates the camera. obdeect side: `sources` → propagation → to-simtel → sim_telarray (camera), one thread. Both sides run the same sim_telarray camera.

| Quantity | Value |
|---|---|
| sim_telarray side: LightEmission wall s, sim_telarray wall s | 0.054, 0.148 |
| obdeect side: wall s | 0.167 |
| **obdeect / sim_telarray side, wall** | **0.83** |
| obdeect / sim_telarray side, wall, fastest repetitions | 0.80 |

| obdeect program | Wall s |
|---|---|
| sources: light source | 0.005 |
| propagation: run (atmosphere) | 0.003 |
| stream: to-simtel (hand the photons to sim_telarray) | 0.003 |
| sim_telarray, camera only | 0.156 |

| Check | Measured | Result |
|---|---|---|
| obdeect wall ≤ LightEmission + sim_telarray wall | 0.83 × sim_telarray | PASS |

All runs completed with their repetitions (2 checks).

<!-- notes:performance-flasher-dish -->
Only about 56 000 photons reach the camera. Both sides are therefore dominated by start-up and by the sim_telarray camera, which takes 0.15 s on each side.
- LightEmission (`ff-1m`) itself takes 0.05 s, mostly start-up.
- The obdeect flasher takes 0.005 s.

LightEmission writes an uncompressed file here (gzip added 0.03 s to its side before).

To compare the light sources themselves, the case needs far more photons per event.
<!-- /notes -->

Data: `performance.flasher-dish` · [metrics.json](performance.flasher-dish/metrics.json) · [details.json](performance.flasher-dish/details.json)

![Wall time per program on both sides.](performance.flasher-dish/wall-by-stage.png)

*Wall time per program on both sides.*

<a id="performance-point-source"></a>

### Distant point source of light — PASS

2 events of 10 million photons each from an isotropic point source 100 km away. sim_telarray side: LightEmission `xyzls` writes an uncompressed photon file (as the obdeect stream), then sim_telarray simulates the camera. obdeect side: `sources` → propagation → to-simtel → sim_telarray (camera), one thread. Both sides run the same sim_telarray camera.

| Quantity | Value |
|---|---|
| sim_telarray side: LightEmission wall s, sim_telarray wall s | 0.045, 0.122 |
| obdeect side: wall s | 0.126 |
| **obdeect / sim_telarray side, wall** | **0.74** |
| obdeect / sim_telarray side, wall, fastest repetitions | 0.75 |

| obdeect program | Wall s |
|---|---|
| sources: light source | 0.002 |
| propagation: run (atmosphere) | 0.001 |
| stream: to-simtel (hand the photons to sim_telarray) | 0.001 |
| sim_telarray, camera only | 0.121 |

| Check | Measured | Result |
|---|---|---|
| obdeect wall ≤ LightEmission + sim_telarray wall | 0.74 × sim_telarray | PASS |

All runs completed with their repetitions (2 checks).

<!-- notes:performance-point-source -->
Same as the flasher. The sim_telarray camera takes 0.12 s on both sides. LightEmission (`xyzls`) takes 0.045 s, the obdeect point source 0.002 s.
<!-- /notes -->

Data: `performance.point-source` · [metrics.json](performance.point-source/metrics.json) · [details.json](performance.point-source/details.json)

![Wall time per program on both sides.](performance.point-source/wall-by-stage.png)

*Wall time per program on both sides.*
