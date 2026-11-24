#!/bin/bash
# Generate a CORSIKA shower input of the performance cases inside the validation image.
#   CORSIKA_DATA=/path/to/corsika7-interaction-tables/interaction-tables \
#     bash obdeect-validation/containers/generate-input.sh gamma-20deg-100-3000GeV-100showers
# Writes obdeect-validation/inputs/<card>.iact and <card>.corsika.log; prints the SHA-256 for
# fixtures/inputs.json. Run from the workspace root, after BUILD_ONLY=1 run-validation.sh.
set -euo pipefail
card=${1:?card name in fixtures/corsika/, without .card}
: "${CORSIKA_DATA:?set CORSIKA_DATA to the CORSIKA interaction-table directory (UrQMD-1.3.1-xs.dat, ...)}"
image=${IMAGE:-obdeect-validation:local}
mkdir -p obdeect-validation/inputs
podman run --rm \
  --volume "$PWD/obdeect-validation:/validation:ro" \
  --volume "$PWD/sim_telarray:/reference:ro" \
  --volume "$PWD/obdeect-validation/inputs:/out:Z" \
  --volume "$CORSIKA_DATA:/corsika-data:ro" \
  --env PYTHONPATH=/validation/src --env CORSIKA_DATA_DIR=/corsika-data \
  "$image" bash -c '
    set -euo pipefail
    card=$1; work=$(mktemp -d); cd "$work"
    python -c "from pathlib import Path; from obdeect_validation.adapters import chain_showers as c; c.prepare_data(Path(\".\"))"
    # CORSIKA reads its tables (qgsdat-II-04, ...) under their original names from the mount.
    sed -e "s|^DATDIR .*|DATDIR /corsika-data/|" -e "s|^TELFIL .*|TELFIL $card.iact|" -e "s|^DIRECT .*|DIRECT ./|" \
        -e "/^EXIT/d" /validation/fixtures/corsika/$card.card > run.card
    echo EXIT >> run.card
    /workdir/simulation_software/corsika7/corsika_qgs2_urqmd_flat < run.card > $card.corsika.log 2>&1
    cp $card.iact $card.corsika.log /out/
  ' _ "$card"
sha256sum "obdeect-validation/inputs/$card.iact" 2>/dev/null || shasum -a 256 "obdeect-validation/inputs/$card.iact"
