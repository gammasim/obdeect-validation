# Validation of the curved and spherical atmosphere

Suite `suites/curvature.json`. It validates the curvature modes of the obdeect atmosphere stage
(`none`, `corsika`, `spherical`; see obdeect `docs/propagation/curvature.md`). Limits are fixed
before results are seen, with a justification in the adapter (`LIMITS`).

## Cases without showers (also in the `propagation` suite)

| Case | Reference | Checks |
|---|---|---|
| `propagation.curvature-chapman` | Chapman integral of an exponential atmosphere, table-form slant integral | transmission of the `corsika` mode against both; flat limit; curved never below plane-parallel at large zenith |
| `propagation.curvature-geometry` | analytic spherical geometry | emission altitude, arrival point and time, dropped photons, flat limit, d²/(2R) shift |
| `propagation.refraction-curved-ray` | independent numpy Bouguer ray tracer (own spline, own quadrature) | shift, direction scale, time and straight-line delay of the `spherical` mode with a profile; flat limit |
| `propagation.curvature-input-guards` | documented mode rules | accepted and rejected combinations through the CLI and the Python API, with the offending key in the message |

The limits for float32 outputs include the float32 half-ulp of the stored value. The Chapman case
computes the zenith cosine from the stored `cx`, `cy` as the stage does.

## Cases on CORSIKA CURVED showers

Inputs: gamma showers at 65 and 75 deg from the curved CORSIKA binary (cards in
`fixtures/corsika/*-curved.card`), and a flat 65 deg control. The cards use the long bunch format
(`TELFIL +`) and full-size telescope spheres. Generate them with
`python -m obdeect_validation.adapters.curved_showers run --geometry curved --zenith 65 ...`;
`run` checks the CORSIKA log for `CURVED VERSION WITH SLIDING PLANAR ATMOSPHERE`, so a flat binary
cannot pass as curved. Hashes are in `fixtures/inputs.json` (files generated 2026-11-24).

| Case | Checks |
|---|---|
| `stream.corsika-curved-flag` | `curvature_applied` of the stream against the independently read EVTH word 79 and bit `0x40` of word 77; refraction bit; log banner |
| `propagation.curvature-corsika-geometry` | slant depth of the median emission altitude in curved geometry within the gamma Xmax window; mean zenith of near-core bunches equals the primary zenith (detector frame); long bunch format; header flags |
| `propagation.curvature-corsika-extinction` | `corsika`-mode transmission against an independent curved computation from the bunches (5 standard errors, numerical floor); the predicted curved-to-plane-parallel ratio is at least 1 and rises with zenith |
| `chain.curvature-high-zenith` | detected p.e. ratio of the obdeect chain to sim_telarray against the predicted ratio; sim_telarray's own transmission against the plane-parallel prediction |

### sim_telarray is plane-parallel

sim_telarray uses the airmass $1/|c_z|$ also for CORSIKA CURVED input (`sim_telarray.c`,
`sim_absorb.c`). The chain case therefore does not expect zero difference: it tests the predicted
ratio of curved to plane-parallel transmission, weighted by quantum efficiency and mirror
reflectivity.

### What the showers can resolve

At 75 deg the curved correction is 2 to 3 % in the airmass, but only about 0.2 % in the mean
transmission, because the photons come from emission heights with a moderate optical depth. A 5
sigma detection of 0.2 % needs about $2.5\times10^7$ bunches, more than the shower files provide.
The chain case and the extinction case are therefore consistency checks within their 1 % floors
and numerical floors; the analytic cases above carry the resolving power.

## Not a case

The comparison with the CORSIKA tables `DISTIP`/`TOFIP` needs an instrumented CORSIKA run; see
`docs/curvature-corsika-tables.md`.

## Known limits of the checks

- The shower cases were written without CORSIKA CURVED output at hand: the energy windows (300 to
  5000 GeV at 65 deg, 1 to 10 TeV at 75 deg) and the telescope altitude (observation level, sphere
  centres at z = 0 as in the flat card) are unverified until the first run.
- A missing shower file ends the case with an error status, not with "blocked".
