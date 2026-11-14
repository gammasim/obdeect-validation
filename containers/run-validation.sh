#!/usr/bin/env bash
set -euo pipefail

SIMTOOLS_PROD_IMAGE=${SIMTOOLS_PROD_IMAGE:-ghcr.io/gammasim/simtools-prod:v0.39.0-v78010-v2025-11-30-rc-generic@sha256:41b5f306407f7cd9f485b4e7da7e0ccb992c4737cd840119c4be3be81cfe9d1b}
if [[ ! "$SIMTOOLS_PROD_IMAGE" =~ @sha256:[0-9a-f]{64}$ ]]; then
  echo "SIMTOOLS_PROD_IMAGE must end with @sha256: and 64 hex digits" >&2
  exit 2
fi
suite=${1:-sources}
output=${2:-obdeect-validation/results/container}
if [[ ! "$suite" =~ ^[a-z0-9-]+$ || ! -f "obdeect-validation/suites/$suite.json" ]]; then
  echo "unknown suite: $suite" >&2
  exit 2
fi
image=obdeect-validation:local

mkdir -p "$output"
output=$(cd "$output" && pwd)
podman build \
  --ignorefile obdeect-validation/containers/context.ignore \
  --build-arg "SIMTOOLS_PROD_IMAGE=$SIMTOOLS_PROD_IMAGE" \
  --file obdeect-validation/containers/Dockerfile.simtools-prod \
  --tag "$image" .
image_id=$(podman image inspect --format '{{.Id}}' "$image")
python3 -c 'import json,sys; json.dump({"base_image":sys.argv[1],"derived_image_id":sys.argv[2]},open(sys.argv[3],"w"),indent=2); open(sys.argv[3],"a").write(chr(10))' \
  "$SIMTOOLS_PROD_IMAGE" "$image_id" "$output/container.json"
# Same layout as the development container: checkouts read-only, results writable.
models=${VALIDATION_MODELS_DIR:-simulation-models-dev}
mounts=(--volume "$PWD/obdeect-validation:/validation:ro"
        --volume "$PWD/sim_telarray:/reference:ro"
        --volume "$(cd "$models" && pwd):/models:ro")
podman run --rm "${mounts[@]}" \
  --volume "$output:/results:Z" \
  "$image" \
  obdeect-validate run \
    --suite "/validation/suites/$suite.json" \
    --environment /validation/environments/workspace.json \
    --output /results --non-interactive
