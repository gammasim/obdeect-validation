# CORSIKA CURVED tables (DISTIP, TOFIP) against obdeect spherical refraction

Status: design note. No case is registered for this comparison yet.

## Decision

The CORSIKA table routine can be built and run standalone without running CORSIKA. The
comparison of obdeect's spherical bending and time with the tables is not yet a valid test: the
arguments of the CORSIKA tables are apparent quantities of the detector's tangent frame, obdeect's
are vacuum spherical quantities, and the sources do not fix the mapping between them. A case would
encode an unverified convention, so none was written. The smallest instrumented CORSIKA run that
settles the mapping is given in section 4.

## 1. What CORSIKA does

`corsika7/src/corsika.F` (CURVED option, `#if __CURVED__`):

- `INRTAB` (line 78229) fills `DISTEF(0:360, 0:452)` and `TOF(0:360, 0:452)` (`RTABLE`
  common, `corsika.h` line 5072). Grid: theta 0 to 90 deg in 0.25 deg steps, height 0 to 113 km in
  0.25 km steps; entries below the observation level or above 89 deg are zero.
- Each entry integrates one ray in 1000 steps: start at the emission height, apparent angle
  `THEAP`, Snell's law `sin(THAPP) = (n_e/n) sin(THEAP)` with `n - 1 = ETADSN * RHOF(h)`,
  `DDIST = DHAPP tan(mean angle)` and `DT = DLEN n / c` (ns).
- `RHOF` (line 30770) is the five-layer density parametrisation of the active atmosphere
  (`HLAY`, `AATM`, `BATM`, `CATM`, `DATM`), set from `BLOCK DATA BLOCK1` (line 6324) and the
  model selection (around line 32880, `MODATM`, `LAYNO`).
- `DISTIP(THEAP, HEAPP)` (line 77969) and `TOFIP` (line 79083) bilinearly interpolate the tables
  (cell index `INT(4 * ...)`, linear in theta and height).
- Use in `CERENK` (lines 77440 to 77600): `RDIST = DISTIP(WEMIS, ZAPP)` and
  `CARTIM = TEMIS * 1e9 + TOFIP(WEMIS, ZAPP)`, with `ZAPP = ZHEM` in the detector frame
  (lines 77360 to 77366). `ZHEM` is the emission height seen from the detector.

## 2. Standalone reference (built and run)

The table routine has no dependency on the rest of CORSIKA beyond its common blocks, so it was
compiled in a scratch directory with the local gfortran (Homebrew GCC 16.2):

1. Copy `corsika.h` next to the extracted units.
2. Extract `BLOCK1` (lines 6324 to 6756), `INRTAB` (78229 to 78369) and `RHOF` (30770 to 30820)
   from `corsika.F` into one file. Each unit carries its own `#define __XXXINC__` lines and
   `#include "corsika.h"`.
3. Replace the three build macros that CORSIKA's build fills in: `__CVERSION__` (7.8010D0),
   `__ICDATE__` (20260101) and `__CDATE__` (a string). These only appear in `DATA` statements.
4. Compile with `gfortran -cpp -ffixed-line-length-none -w -c`.
5. A driver (its own `PROGRAM`, same `#define` set and `#include "corsika.h"`) sets the model
   selection as the CORSIKA routine does for `MODATM = 1`: `HLAY(L) = HLAY0(L, LAYNO(1))`,
   `AATM`, `BATM`, `CATM` from model 1, `DATM = 1/CATM`; then `ETADSN = 0.283D-3 CATM0(1,1) /
   BATM0(1,1)`, `OBSLEV(1) = 215600` (cm), calls `INRTAB` and writes `DISTEF` and `TOF`.

Result: 163,533 rows (361 by 453). Checks:

- `C(1) = 637131500` cm (Earth radius), `C(25) = 2.99792458e10` cm/s.
- Vertical rays: `DISTEF(0, M) = 0` at every height.
- Vertical from 10 km to the observation level: `TOF(0, 40) = 26168.8054 ns`, that is a vacuum
  flight time of 26164.77 ns plus an air delay of 4.04 ns. An independent estimate with the
  atmprof1 density and the same `ETADSN` gives 3.96 ns. The 2 % difference is the five-layer
  density: `RHOF(0) = 1.2298e-3 g/cm3` against 1.1668e-3 in atmprof1 (+5.4 %), and at 10 km
  4.2614e-4 against 4.2066e-4 (+1.3 %).
- Sample values: `DISTEF(340, 40) = 8.8930e6 cm` (85 deg, 10 km), `TOF(340, 40) = 2.9783e5 ns`.

These checks establish that the built routine reproduces the physics of its own atmosphere. They
do not establish the convention needed to compare with obdeect.

## 3. Why the comparison is not yet valid

1. Argument meaning. `DISTEF` is labelled with apparent quantities (`APPARENT EMISSION ANGLE`,
   `APPARENT EMISSION HEIGHT`). It builds points with `DHELP = (HEAPP - OBS) tan(THEAP)` and
   `HE = sqrt(DHELP^2 + (R + HEAPP)^2) - R`, that is a straight line in the detector's
   tangent frame, not a straight line in the Earth-centred sphere that obdeect uses (vacuum
   impact parameter `b = r_t sin(theta_s)`). The two straight lines differ at order distance/R:
   at 85 deg and 10 km the horizontal distance is 0.7 to 0.9 km, so the angle difference is
   about 1.4e-4 rad (0.008 deg), the order of the curvature effect under test.
2. Angle reference. `THEAP` is a cosine of an emission angle in the detector frame (`WEMIS` after
   the Earth-angle rotation). obdeect's input is the vacuum zenith at the telescope, `theta_s`.
   A candidate mapping is `THEAP = theta_s` with `HEAPP` chosen so that the true altitude
   `HE(HEAPP, theta_s)` equals obdeect's emission altitude `h_e`. This is a hypothesis, not a
   derivation from the sources, and it must be checked before use.
3. Shift and time conventions. `DISTEF` is the horizontal distance travelled by the ray (accumulated
   as `DDIST`, positive in the table). obdeect's `shift` is zero or negative along the direction
   of motion. `TOFIP` is the full flight time from the apparent emission point (`CARTIM` adds
   `TEMIS`, the time to that point), while obdeect's `delta_t` is the excess over the vacuum
   straight line. The comparison must use the same excess: `TOF(theta, H) - straight-line time`.
4. Atmosphere. obdeect's refraction uses the profile in its input. For a like-for-like comparison
   the profile must be written from `RHOF` (density) and `ETADSN * RHOF` (index minus one), so
   that both sides use the five-layer model. With atmprof1 the index differs by about 4.5 % at
   sea level; the effect on the shift is not quantified here.
5. Grid error. The table is linear in theta (0.25 deg) and height (0.25 km). Its interpolation
   error is not yet measured. The limit must come from that measurement, made before the
   comparison: evaluate the same integral directly (the `INRTAB` inner loop at finer steps) at
   cell midpoints (0.125 deg, 0.125 km) and compare with the bilinear value. Do not set the limit
   from the observed agreement.

## 4. Smallest instrumented CORSIKA run

Purpose: settle items 1 to 3 with the code that actually runs, without changing the physics.

- Input: CURVED build (`corsika7-config/config/78010/config_*_curved.h` or the image's curved
  variant), one shower (`NSHOW 1`), one primary gamma at 85 deg and 1 TeV, one telescope at the
  core (sphere radius from `curved_showers.layout(85)`), observation level 215600 cm, atmosphere 1.
- Enable the existing CERENK debug output with the `CDEBUG T` card (`corsika.F` line 12059).
  It writes `CERENK: WEMIS,ZEM,STHE,RDIST=` at line 77457 (`RDIST = DISTIP(WEMIS, ZAPP)`). It prints
  `ZEM` but not `ZAPP`, so a patch is needed: add `ZAPP` to that write and print `TOFIP` next to
  `CARTIM` (lines 77508, 77523, 77588, 77600). The patch affects output only.
- Collect about 100 photons. For each: `WEMIS`, `ZAPP`, `RDIST`, `TOFIP`.
- Compare with the standalone routine at the same `(WEMIS, ZAPP)`: `DISTIP` and `TOFIP` must
  agree to interpolation precision (the same grid), which checks the build, not the physics.
- Only then compare with obdeect at the mapping of item 2, using the excess time of item 3.

## 5. Files and provenance

- The standalone harness was built in the analysis scratch area, not in this repository. The
  recipe above reproduces it: the line ranges are those of the CORSIKA 7 source tree in
  `simpipe/simulation_software/corsika7/src/corsika.F`.
- Nothing in this note changes obdeect or the validation limits.
