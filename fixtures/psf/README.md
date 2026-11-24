# PSF fixture: PSFcurve_data_v2.ecsv

Provenance of `PSFcurve_data_v2.ecsv`, used by `optics-cumulative-psf-angles`
(`fixtures/optics-inputs/optics-cumulative-psf-angles.json`, adapter `optics_psf_angles`).

- **What it is:** measured point-spread function of LSTN-01 at zenith about 20 deg, with Vega
  (ECSV header comment `PSF taken at Zd = ~20[deg], with Vega`). Columns: `radius` (mm),
  `differential_value` and `integral_value` of the PSF. It is a simulation-model parameter
  (`setting_psf_parameters` of LSTN-01, description "PSF measurement as function of offsets angle").
- **Source repository:** `git@github.com:gammasim/simtools-tests.git` (simtools integration-test
  resources), checked out at `/Users/maierg/DESYCloud/CTA/Pipelines/MC-Pipeline/gammasim-tools/simtools-tests`.
- **Path in that repository:** `simtools-tests/v0.39.0/integration_tests/downloaded/PSFcurve_data_v2.ecsv`
  (the v0.39.0 test release; it was downloaded there from the model-parameter repository path
  `input/LSTN-01/setting_psf_parameters/019d7772-b58f-78d4-9d7e-cd7fcbc2ccc9/PSFcurve_data_v2.ecsv`
  listed in `v0.39.0/integration_tests/config_files/download_files.yml`).
- **Version and commit:** tracked in `simtools-tests` at commit `11c35202625fd69d39c4300148ba82360093722d`
  ("v0.39.0 initial release of tests"), blob `8bc20b49771d0103db0372c0fdf250209fcc5e25`. The same file
  (identical SHA-256) is also present in the v0.34.0 to v0.40.0 integration tests and in `legacy/v0.32.0`.
- **SHA-256:** `aedc136b77da6fbcf3ca789908fddc408204f8a4e24a2b0e15eead4d79294e33`
  (check with `shasum -a 256 PSFcurve_data_v2.ecsv`).

The file is an unmodified copy; obdeect-validation does not change it.
