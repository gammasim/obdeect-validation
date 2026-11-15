# CORSIKA CURVED tables (DISTIP, TOFIP) against obdeect spherical refraction

Status: the analysis and the output patch are prepared and tested without CORSIKA. The
instrumented run has not been made: it needs the owner's CORSIKA build. The case
`propagation.curvature-corsika-tables` (`cases/propagation/curvature-corsika-tables.json`) is not in
a suite yet.

## Decision

The CORSIKA table routine can be built and run standalone without running CORSIKA, and the CURVED
build can be traced with output only (section 4). The comparison of obdeect's spherical bending and
time with the tables is still not a valid test. The arguments of the CORSIKA tables are apparent
quantities of the detector's frame, obdeect's are vacuum spherical quantities, and the sources do
not fix the mapping between them. The analysis therefore reports residuals under an explicit
hypothesis, without limits. The first real trace decides whether the hypothesis holds.

## 1. What CORSIKA does (checked against `corsika7/src/corsika.F`)

Line numbers refer to the checkout at commit 4363a29c (blob 4548f67a, unmodified there).

- `INRTAB` (line 78229) fills `DISTEF(0:360, 0:452)` and `TOF(0:360, 0:452)` (`RTABLE` common).
  Grid: theta 0 to 90 deg in 0.25 deg steps, height 0 to 113 km in 0.25 km steps. Entries at or
  below the observation level, or above 89 deg, are zero. Each entry integrates one ray in 1000
  steps from the apparent emission point down to the observation level: Snell's law
  `sin(THAPP) = (RHE/RHST) sin(THEAP)`, `DDIST = DHAPP tan(mean angle)` (horizontal distance) and
  `DT = DLEN n / c` (ns).
- Refractive index in `INRTAB` (line 78302): `IATMOX >= 1` uses `REFIDX` (the external atmosphere
  of the bernlohr code, set by the ATMOSPHERE card's first argument); otherwise `1 + ETADSN RHOF(h)`
  with the five-layer density.
- `DISTIP(THEAP, HEAPP)` (line 77969) and `TOFIP(THEAP, HEAPP)` (line 79083) interpolate
  bilinearly. `THEAP2 = ACOS(THEAP) FACT 4` (theta in units of 0.25 deg), `HEAPP2 = HEAPP 1e-5 4`
  (height in units of 0.25 km), `N2 = INT(THEAP2)`, `M2 = INT(HEAPP2)`. Below the observation level
  the weights use `HM2 = OBS` (`IF (HM2 .LT. OBS)` at lines 78029 and 79146). `TOFIP` writes a line
  to `MDEBUG` on every call when `LCERDB` is set (line 79160).
- `CERENK` (line 76850), CURVED branch (`__CURVED__`):
  - `ZAPP = ZHEM` when `DETSYS` is false, `ZAPP = ZEM` otherwise (lines 77361-77366).
  - `RDIST = DISTIP(WEMIS, ZAPP)` (line 77456) is called unconditionally in the CURVED branch. The
    photon's ground point is `XCER2 = XEMIS + RDIST/STHE * UEMIS2` (line 77461): `RDIST` is the
    full horizontal offset from the emission point's vertical, not a shift relative to a straight
    line. In the flat branch the same offset is `PATHCR * STHE = (ZEM - OBS) tan(theta)`.
  - `CARTIM = TEMIS*1e9 + TOFIP(WEMIS, ZAPP)` (lines 77508 and 77523 in the IACT branch; 77588,
    77600, 77664 and 77676 in the non-IACT branch). Both arms of every `IF (FREFRX)` use `TOFIP`.
  - The `CDEBUG` output (line 77457) prints `WEMIS, ZEM, STHE, RDIST` but not `ZAPP`, and nothing
    next to `CARTIM`.
- Refraction flag. `ATMOSPHERE IATM Y` sets `IATMOX = IATM` and `FREFRX = .TRUE.` (lines 12005-12006;
  default `FREFRX = .FALSE.`, line 11601). `FREFRX` sets EVTH word 77 bit 0x10 (line 33806) and
  selects the FREFRX arm of CERENK. In a CURVED build with `__ATMEXT__` both arms compute the same
  `CARTIM` and `RDIST` from the tables. So the tables are used with and without `FREFRX`; the
  refraction flag does not switch them on. The earlier note said the tables are used only with
  ATMEXT refraction; the source does not show that. The first argument of the card (`IATMOX`)
  chooses the index model in `INRTAB`.

## 2. Standalone reference (recorded earlier, not re-run for this note)

The table routine has no dependency on the rest of CORSIKA beyond its common blocks, so it was
compiled in a scratch directory with the local gfortran (Homebrew GCC 16.2):

1. Copy `corsika.h` next to the extracted units.
2. Extract `BLOCK1` (lines 6324 to 6756), `INRTAB` (78229 to 78369) and `RHOF` (30770 to 30820)
   from `corsika.F` into one file. Each unit carries its own `#define __XXXINC__` lines and
   `#include "corsika.h"`.
3. Replace the three build macros that CORSIKA's build fills in: `__CVERSION__` (7.8010D0),
   `__ICDATE__` (20260101) and `__CDATE__` (a string). These only appear in `DATA` statements.
4. Compile with `gfortran -cpp -ffixed-line-length-none -w -c`.
5. A driver sets the model selection as the CORSIKA routine does for `MODATM = 1`: `HLAY(L) =
   HLAY0(L, LAYNO(1))`, `AATM`, `BATM`, `CATM` from model 1, `DATM = 1/CATM`; then `ETADSN =
   0.283D-3 CATM0(1,1) / BATM0(1,1)`, `OBSLEV(1) = 215600` (cm), calls `INRTAB` and writes `DISTEF`
   and `TOF`.

Recorded result: 163,533 rows (361 by 453). `C(1) = 637131500` cm, `C(25) = 2.99792458e10` cm/s.
Vertical rays give `DISTEF(0, M) = 0`. Vertical from 10 km to the observation level: `TOF(0, 40) =
26168.8054 ns` (vacuum 26164.77 ns plus 4.04 ns air delay; an atmprof1 estimate gives 3.96 ns,
the 2 % difference is the five-layer density). Sample: `DISTEF(340, 40) = 8.8930e6 cm`,
`TOF(340, 40) = 2.9783e5 ns`. These checks show that the routine reproduces the physics of its own
atmosphere. They do not establish the convention needed to compare with obdeect.

## 3. Why the comparison is not yet valid

1. Argument meaning. `DISTEF` is labelled with apparent quantities. It builds points with
   `DHELP = (HEAPP - OBS) tan(THEAP)` and `HE = sqrt(DHELP^2 + (R + HEAPP)^2) - R`: the horizontal
   offset is taken in the tangent frame, the altitude on the sphere. obdeect's geometry is the
   vacuum sphere (impact parameter `b = r_t sin(theta_s)`). The two differ at order distance/R.
   At 85 deg and 10 km the difference is about 1.4e-4 rad.
2. Angle reference (hypothesis H, not verified). `THEAP` is a cosine of an emission angle in the
   detector frame (`WEMIS` after the Earth-angle rotation). obdeect's input is the vacuum zenith at
   the telescope, `theta_s`. H: `THEAP = theta_s`, and the apparent height `HEAPP` satisfies
   `HE(HEAPP, theta_s) = h_e`. Nothing in the sources fixes H.
3. Which emission altitude is true. `ZEM` is the altitude of the step in CORSIKA's intermediate
   frame after the `CDDIF` transform (lines 77116-77118: `ZBEG = (ZBEG + R)/CDDIF - R`). The
   altitude before that transform is not printed. H needs `h_e` from `ZEM`, so this is open.
4. The horizontal offset. `RDIST` contains the straight part `(h_e - h_obs) tan(theta)` as well as
   the refraction shift (section 1). The analysis therefore compares `RDIST` with
   `(h_e - h_obs) tan(theta_s) + shift`, with the signed shift from obdeect. An earlier version
   compared `RDIST` with `|shift|` alone; that is dimensionally inconsistent (the difference is
   about 20 km at 70 deg and 8 km).
5. Shift and time conventions. `TOFIP` is the flight time of the refracted apparent path from the
   emission point to the observation level. The comparable obdeect quantity is `t = s/c + excess`
   with `s = (h_e - h_obs)/cos(theta_s)` (the tangent-plane straight path, as in obdeect's
   reference). So `excess_CORSIKA = TOFIP - s/c`.
6. Atmosphere. `IATMOX = 1` uses `REFIDX` (the bernlohr external atmosphere). The obdeect reference
   uses the atmprof table of the sim_telarray configuration (`ATMOSPHERE_PROFILE`). The index
   models differ; the effect on the shift is not quantified here (an earlier estimate of the index
   difference at sea level is 4.5 %).
7. Grid error. The table is linear in theta (0.25 deg) and height (0.25 km). The analysis checks the
   interpolation arithmetic only (section 5). The error of the grid against the exact integral is
   not measured: evaluate the `INRTAB` integral at cell midpoints (0.125 deg, 0.125 km) with finer
   steps and compare. The limit must come from that measurement, not from the observed agreement.

## 4. Instrumented run: exact steps

### 4.1 Patch the CORSIKA source

The patch is `fixtures/corsika-patches/curvature-tables-debug.patch` in obdeect-validation. It
changes output only: a new function `CVTBUN` (the trace unit), a new subroutine `CVTBHD` called at
the end of `INRTAB`, and a trace block in `CERENK` after the IACT `TELOUT` call. It adds no card
keyword, no common block and no physics change. Its SHA-256 starts `2ba42d8149bbd1d1`.

The patch is a unified diff whose target is `corsika.F`, so it is applied with `-p0` from the
directory that holds `corsika.F` (`src/` of the corsika7 checkout):

```sh
cd <corsika7>/src                                   # the source of the CURVED build
shasum -a 256 corsika.F                             # expect 62b0859db3cadc60e137... (4363a29c)
patch -p0 --dry-run < <validation>/fixtures/corsika-patches/curvature-tables-debug.patch
patch -p0 < <validation>/fixtures/corsika-patches/curvature-tables-debug.patch
shasum -a 256 corsika.F                             # expect 6c0568edb8266c2e00f6... after the patch
```

`git apply --check` also accepts the patch, and a reverse dry run succeeds. The checkout's
`src/corsika.F` reports `CORSIKA 7.8041`, while `corsika7-config` has headers for 78010 and 78050.
If the build's source is another version, the patch may apply with offsets or fail; check with
`--dry-run` first and do not force it.

Front-end check, done here: `gfortran -cpp -DHAVE_CONFIG_H -fsyntax-only` on the original and on
the patched file with the curved configuration (section 4.2) passes for both, with the same warning
count (10, the `H` format specifier). This is a parse check, not a build; no binary was made.

### 4.2 Build the curved binary

- Configuration: `simpipe/simulation_software/corsika7-config/config/78010/config_qgs2_urqmd_curved.h`
  (QGSJetII, UrQMD; `__CURVED__ 1`, `__CERENKOV__ 1`, `__IACT__ 1`, `__ATMEXT__ 1`). Other choices
  in the same directory are possible; the production choice is the owner's.
- The source includes `config.h` when `HAVE_CONFIG_H` is defined (`corsika.F` line 197), so the
  build must provide `config.h` (the coconut step copies the header) and compile with
  `-DHAVE_CONFIG_H`. Some macros in the header have no value (for example `#define __NRREXT__`), and
  `#if __NRREXT__` then does not parse; the front-end check gave such macros the value 1. How the
  production build handles this is not visible in this checkout.
- The image's build command is in simtools (`docker/Dockerfile-corsika7`), not in this checkout.
  Use the same command as for the production curved binary, with the patched source.
- Record the binary path, its SHA-256 and the SHA-256 of the patched `corsika.F`.

### 4.3 Card

One shower, one 1 TeV gamma, one telescope at the core, long bunch format (`TELFIL +`), sphere
radius 1200 cm (the project's decision for the long format; the compact-format radius limit of
`curved_showers.py` does not apply to `+`). Written in the style of
`fixtures/corsika/gamma-65deg-300-5000GeV-curved.card`; the 75 deg card differs only in `THETAP`
and the `TELFIL` name.

```text
RUNNR 1
EVTNR 1
NSHOW 1
PRMPAR 1
ERANGE 1000 1000
ESLOPE -2.0
THETAP 65 65
PHIP 0.0 0.0
SEED 17 0 0
SEED 29 0 0
SEED 43 0 0
SEED 61 0 0
OBSLEV 215600.0
ATMOSPHERE 1 Y
MAGNET 20.0 40.0
ARRANG 0.0
HADFLG 0 0 0 0 0 2
ECUTS 0.3 0.3 0.02 0.02
ELMFLG F F
CSCAT 1 0.0 0.0
DATDIR <corsika data>/
DIRECT ./
TELESCOPE 0.0 0.0 1200.0 1200.0
TELFIL +./gamma-65deg-1TeV-curved.iact
IACT SPLIT_AUTO 0
IACT IO_BUFFER 800MB
IACT MAX_BUNCHES 1000000
EXIT
```

Points to check in the output, not assumed here:

- `ERANGE 1000 1000` gives a fixed energy: `LLIMIT .EQ. ULIMIT` sets `ISPEC = 0` (line 18680,
  "PRIMARY ENERGY IS FIXED AT ... GEV").
- `ATMOSPHERE 1 Y`: `IATMOX = 1` (index from `REFIDX`) and `FREFRX = T` (EVTH word 77 gets 16).
- `CSCAT 1 0.0 0.0` is not verified here (the fixtures use `CSCAT 5 30000.0 0.0` for three
  telescopes). Check that the impact point of the shower in the IACT file is at the telescope.
- `DATDIR`: the CORSIKA data directory (the environment file sets `CORSIKA_DATA_DIR` to
  `/workdir/simulation_software/corsika7`).

### 4.4 Run

Use an empty run directory. Do not put `CDEBUG` or `DEBUG` in the card: the existing `LCERDB`
prints (for example in `TOFIP`, on every call) would flood the log.

```sh
mkdir run-65 && cd run-65
export OBDEECT_CURVTAB=1                       # the trace is written to ./curvature-tables.dat
<corsika binary> < ../gamma-65deg.card > corsika.log 2>&1
grep "CURVED VERSION WITH SLIDING PLANAR ATMOSPHERE" corsika.log   # curved binary
grep "PRIMARY ENERGY IS FIXED" corsika.log                         # 1000 GeV
head -1 curvature-tables.dat                   # CURVTAB 1 ...
grep "^CONST" curvature-tables.dat             # C1_CM, C25_CMPS, OBSLEV_CM, ETADSN, MODATM, IATMOX, FREFRX
grep -c "^ROW" curvature-tables.dat            # 32
grep -c "^BUNCH" curvature-tables.dat          # > 0
```

The trace file lists at most 200 bunches per shower (`NCVTMX`), those that hit a telescope
(`ITHIT .NE. 0`), in the order CORSIKA produces them. Each `BUNCH` line holds the shower number,
the bunch number, `N2`, `M2` and 22 reals in the order of the patch format (WEMIS, ZAPP, ZEM, STHE,
RDIST, TOFIP, TEMIS*1e9, CARTIM, UMEAN, VMEAN, WMEAN, UEMIS2, VEMIS2, CTHETA, four DISTEF and four
TOF corners of the cell). `ROW` lines are `N = 4 theta(deg)`, `M = 4 h(km)`, `DISTEF` and `TOF`
for theta 0, 30, 60, 70, 75, 80, 85, 89 deg and h 3, 8, 20, 60 km.

### 4.5 Collect

For each zenith, keep:

- `curvature-tables.dat` (renamed to `inputs/curvature-tables-<zen>deg.dat`, for 65 and 75 deg),
- `corsika.log`, the card, the IACT file, the binary and the patched `corsika.F` hashes,
- the commands and the environment (`OBDEECT_CURVTAB` set, `CDEBUG` not set).

Keep these outside git.

### 4.6 Analyse

```sh
VALIDATION_INPUTS=<dir with the two trace files> ATMOSPHERE_PROFILE=<atmprof1.dat> \
  python -m obdeect_validation.adapters.curvature_corsika_tables
```

The case is run in the usual way (`obdeect-validate run --case cases/propagation/curvature-corsika-tables.json
--environment environments/workspace.json ...`). Missing trace files end with exit code 3. The
metrics are the interpolation consistency (section 5); the mapping residuals are in `details.json`
and in `mapping_residuals.png`, without limits.

## 5. What the analysis checks, and what it does not

Checked (metrics, exact limits):

- `table_interpolation_max_relative_error <= 1e-6`. DISTIP and TOFIP recomputed from the printed
  cell corners with the Fortran arithmetic (INT indices, the height clamp to the observation level,
  the two-step interpolation). Relative to the largest corner of the cell; the printed values have 15
  digits, so the agreement should be near 1e-13.
- `cell_index_mismatch_count == 0`. The indices recomputed from the printed WEMIS and ZAPP equal the
  printed `N2`, `M2`.
- `bunch_count_positive == 1` and `file_count == 2`.

Also in details, not gated: `CARTIM = TEMIS*1e9 + TOFIP` (`cartim_max_relative_error`), and the
printed grid rows against the cell corners where a corner is a printed node.

Not checked:

- The table against the exact integral (item 7 of section 3).
- The mapping (section 3, items 1 to 6). It is reported as residuals, in bins of zenith and emission
  height, under hypothesis H, with no limits.

## 6. Open points for the first real run

1. Confirm H: `theta_s = arccos(WEMIS)` and `h_e = ZEM`. The altitude residual
   `HE(ZAPP, theta) - ZEM` tests the height part; zero means the apparent geometry of INRTAB matches
   `ZEM`. If it is not zero, decide whether `ZEM` or the pre-transform altitude is the true one (section
   3, item 3). The trace does not print the pre-transform altitude; add it to the patch if needed.
2. Confirm the meaning of `RDIST` (section 3, item 4) from the dependence of the shift residual on zenith
   and height.
3. Put both models on the same index. `IATMOX = 1` uses `REFIDX`; the obdeect side uses atmprof1.
4. Measure the grid error (section 3, item 7) before setting any limit on the mapping.
5. Bunch sample. The first 200 bunches per shower are the earliest on a telescope. Decide whether to
   sample across the shower (a change to `NCVTMX` or to the selection in the patch).
6. Core position. Check the impact point in the IACT file (`CSCAT 1 0.0 0.0` is not verified).
7. The obdeect side is the numpy reference (`refraction_curved_reference.ray_integrals`), not the
   obdeect `propagate` API. The existing stream helper places the telescope at 2000 m, the CORSIKA
   runs at 2156 m. A public-API path needs a site parameter. The reference is the one the
   refraction case already validates against obdeect.
8. CORSIKA version: the checkout reports 7.8041; the configuration headers are for 78010 and 78050.
   Check that the patch applies to the version of the production binary.
9. Inputs. `cases/` cannot list inputs (the schema has no `inputs` key). Add the two trace files to
   the `inputs` list of the environment file (as the shower files are) so that a missing file shows as
   BLOCKED. This is the owner's decision; the suites were not changed here.
10. Add the case to a suite after the first run.

## 7. Files and provenance

- Patch: `fixtures/corsika-patches/curvature-tables-debug.patch` (against `corsika7/src/corsika.F`
  at 4363a29c). Verified: `patch -p0 --dry-run`, `git apply --check`, a reverse dry run, and a
  reproduction of the expected patched file. A pytest checks the patch against a copy of the
  source when `CORSIKA_SOURCE` is set (skipped otherwise).
- Adapter: `src/obdeect_validation/adapters/curvature_corsika_tables.py`.
- Case: `cases/propagation/curvature-corsika-tables.json`.
- Tests: `tests/test_curvature_corsika_tables.py` (synthetic traces in the patch format; the bilinear
  arithmetic against a hand-computed example, the clamp and a brute-force reference; residual
  binning; the vacuum limit of the mapping).
- Nothing in this note changes obdeect or the validation limits.
