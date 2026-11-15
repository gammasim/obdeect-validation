# obdeect versus sim_telarray, 2026-11-14

Run: `bash obdeect-validation/containers/run-validation.sh full obdeect-validation/results/full-7.0.0`
(suite `full`, 19 cases) in the simtools-prod based container; obdeect `0.8.1.dev52+g6721ba04b`,
sim_telarray v2025-11-30-rc, simulation models 7.0.0 for all cases (optics and chain). Wall time
908 s. Each case directory holds `metrics.json`, `details.json` (limits, parameter mapping,
reference findings) and its figure. Signed values are obdeect / reference − 1. Supersedes the
earlier 2026-11-14 run that used models 6.0.2 for optics.

**Result: 19 PASS, 0 FAIL as run** (`summary.json` status PASS, no interrupted cases). The night-sky
camera ratio against nsbls (0.80) is recorded only, not a gated metric, so it no longer fails the
case (see below).

| Stage | Case | Agreement | Verdict |
|---|---|---|---|
| Stream | CORSIKA → stream → sim_telarray | 300 telescope events: p.e. +0.07 ± 0.09 %, 0 trigger mismatches | PASS |
| Stream | near-threshold showers (20–100 GeV) | 97 events, 291 telescope samples, 0 trigger mismatches | PASS |
| Stream | LED stream → EventIO, strict reader | block order as CORSIKA; 3/3 events, 60/60 photons | PASS |
| Propagation | analytic Beer–Lambert, 0–60° | max transmission error 0 | PASS |
| Propagation | sim_telarray extinction tables | max error 4.5e-8 (tol. 1e-6) | PASS |
| Propagation | refraction/timing vs compiled CORSIKA atmo | displacement 0.44 mm, time 1.6 ps | PASS |
| Chain | obdeect atmosphere + sim_telarray camera | 300 telescope events: p.e. +0.12 ± 0.15 %, image χ² 1.12 | PASS |
| Optics | LST, 0–4° (7.0.0) | D80 ≤ 0.53 %, area ≤ 0.09 %, focal length ≤ 0.02 % | PASS |
| Optics | MST, 0–4° (7.0.0) | D80 ≤ 1.32 %, area ≤ 0.18 %, focal length ≤ 0.06 % | PASS |
| Optics | SST, 0–4° (7.0.0) | D80 ≤ 0.13 %, area ≤ 0.33 %, focal length ≤ 0.01 % | PASS |
| Chain | obdeect optics + sim_telarray camera (Bypass_Optics=1), LSTN-01 7.0.0 | 100 events: p.e. +0.06 ± 0.19 %, image χ² 0.99, time +0.02 ns | PASS |
| Sources | LED, beam scan, laser geometry | exact / photon ratio 1.000 ± 0.006 / path 131.60 vs 131.65 m | PASS |
| Sources | flasher dish vs ff-1m | photons 1.0015 ± 0.0008; camera p.e. 1.0013 ± 0.0019, χ² 1.03 | PASS |
| Sources | point source vs xyzls | accepted photons 0.995 ± 0.012; time +0.004 ns; camera p.e. 1.010 ± 0.010 | PASS |
| Sources | flasher secondary vs ff-gct | photons 1.0002 ± 0.0022; time +0.008 ns; camera p.e. 1.001 ± 0.006, χ² 1.08 | PASS |
| Sources | muon track vs fake-muon and analytic | yield 1.0003 ± 0.003; time vs analytic 0.000 ns, KS 0.0002; camera p.e. 0.997 ± 0.009 | PASS |
| Sources | night sky vs analytic and nsbls | rate exact; area-uniform entry KS 0.0007; flat pulse KS 0.0006; camera p.e. 0.80 ± 0.007 vs nsbls (recorded only) | PASS |

Change against the 6.0.2 run (same cases, same limits):

- Optics: numbers differ per telescope; all within the provisional limits (LST D80 0.53 %, MST
  1.32 %, SST area 0.33 % vs 0.39 % before).
- Stream, chain (atmosphere and optics) and propagation: unchanged within the printed precision,
  except the sim_telarray table error (0 → 4.5e-8, still well inside tolerance).
- Sources: camera p.e. ratios moved by about 0.5–1 % (flasher dish 1.0019 → 1.0013; point source
  0.999 → 1.010; flasher secondary 0.996 → 1.001; muon track 0.981 → 0.997). Night sky unchanged
  (0.80). The cause was not isolated: obdeect moved from `0.8.1.dev31` to `dev52` as well as the model.

Reference findings (recorded, not judged; details in each `details.json`):

- fake-muon: arrival-time distribution differs from the analytic one by KS 0.021 (mean −5 ps).
- nsbls: entry radius `1 − sqrt(u)` instead of `sqrt(u)` (flux not isotropic, too much light
  at the dish centre), and the entry point is rotated with the azimuth of `d` instead of `−d`
  (`IactLightEmission.cc` ~l.1395), which breaks impact parameters and times. The night-sky camera
  ratio (0.80) is therefore recorded only; with the same sampler in both programs (2026-11-13 run)
  it was 0.998 ± 0.009.
- ff-gct: near-horizontal photons (|cz| < 0.3) carry a refraction time term of up to 4 ns.
- sim_telarray skips an event block larger than its input buffer (default limits) and still
  exits 0; the flasher-secondary camera stage therefore uses 1e6 photons per flasher. Not an
  obdeect format problem: with `CORSIKA_IO_BUFFER=2GiB MAX_BUNCHES=50000000` the 1.28 GB
  stream is read completely (568833 p.e., 10x the 1e6 run).
- Point source: raw photon ratio 1.039 ± 0.001 vs xyzls is recorded but not gated; the gated
  quantity is the accepted-photon ratio (0.995).

Limits are set before the runs and recorded in each `details.json`; optics and source limits
are provisional and need scientific approval.

Not covered: `Bypass_Optics=3` (absent), dual-mirror shower chain, scattered laser light,
zenith angles other than 20° in the shower chain, per-mirror and incident-angle optics cases.

![optics chain](chain.optics-bypass1/camera-comparison.png)
![muon track](sources.muon-track/muon-track.png)
![night sky](sources.night-sky/night-sky.png)
![flasher secondary](sources.flasher-secondary/stage1.png)
