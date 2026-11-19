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
# MODE3=1 runs the suite with sim_telarray Bypass_Optics=3 (containers/Dockerfile.mode3, tag obdeect-validation:mode3).
run_image=$image
extra_env=()
# Environment file of the run (environments/<name>.json): workspace (all suites), performance (timing).
environment=${ENVIRONMENT:-workspace}
if [[ ! "$environment" =~ ^[a-z0-9-]+$ || ! -f "obdeect-validation/environments/$environment.json" ]]; then
  echo "unknown environment: $environment" >&2
  exit 2
fi

# One-line description of the host CPU (the container cannot see it): model, cores, podman VM CPUs.
host_description() {
  local text model cores perf eff vm
  if [[ "$(uname)" == Darwin ]]; then
    text=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo "unknown CPU")
    if perf=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null) && eff=$(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null); then
      text="$text ($perf performance + $eff efficiency cores)"
    fi
    vm=$(podman machine inspect --format '{{.Resources.CPUs}}' 2>/dev/null || true)
    if [[ -n "$vm" ]]; then
      text="$text, podman VM $vm CPUs"
    fi
  else
    model=$(lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -n 1 || true)
    cores=$(nproc 2>/dev/null || echo "?")
    text="${model:-unknown CPU}, $cores cores"
  fi
  printf '%s\n' "$text"
}

mkdir -p "$output"
output=$(cd "$output" && pwd)
# SKIP_BUILD=1 reuses the images built by an earlier call (several cases one after another).
if [[ "${SKIP_BUILD:-0}" != 1 ]]; then
  podman build --rm \
    --ignorefile obdeect-validation/containers/context.ignore \
    --build-arg "SIMTOOLS_PROD_IMAGE=$SIMTOOLS_PROD_IMAGE" \
    --file obdeect-validation/containers/Dockerfile.simtools-prod \
    --tag "$image" .
fi
if [[ "${MODE3:-0}" == 1 && "${SKIP_BUILD:-0}" != 1 ]]; then
  # sim_telarray source is the git tree at SIMTEL_REF (default: the mode-3 commit), exported with git archive.
  simtel_ref=${SIMTEL_REF:-b8f2494}
  simtel_commit=$(git -C sim_telarray rev-parse --verify "$simtel_ref^{commit}")
  simtel_src=$(mktemp -d "${TMPDIR:-/tmp}/simtel-src.XXXXXX")
  git -C sim_telarray archive "$simtel_commit" | tar -x -C "$simtel_src"
  # hessio at HESSIO_REV (pinned), exported with git archive from the clone HESSIO_REPO.
  hessio_rev=${HESSIO_REV:-604578331f4696703001e4af3b95a19cdcaa3852}
  : "${HESSIO_REPO:?set HESSIO_REPO to a hessio git clone containing $hessio_rev}"
  hessio_src=$(mktemp -d "${TMPDIR:-/tmp}/hessio-src.XXXXXX")
  git -C "$HESSIO_REPO" archive "$hessio_rev" | tar -x -C "$hessio_src"
  podman build --rm \
    --ignorefile obdeect-validation/containers/context.ignore \
    --build-arg "VALIDATION_IMAGE=$image" \
    --build-arg "SIMTEL_REF=$simtel_commit" \
    --build-arg "HESSIO_REV=$hessio_rev" \
    --build-context "simtel=$simtel_src" \
    --build-context "hessio=$hessio_src" \
    --file obdeect-validation/containers/Dockerfile.mode3 \
    --tag obdeect-validation:mode3 .
fi
if [[ "${BUILD_ONLY:-0}" == 1 ]]; then
  exit 0
fi
if [[ "${MODE3:-0}" == 1 ]]; then
  run_image=obdeect-validation:mode3
  extra_env=(--env SIMTEL_CANDIDATE_BINARY=/opt/mode3/sim_telarray)
fi
image_id=$(podman image inspect --format '{{.Id}}' "$run_image")
host=$(host_description)
extra_env+=(--env "OBDEECT_VALIDATION_HOST=$host")
python3 -c 'import json,sys; json.dump({"base_image":sys.argv[1],"derived_image_id":sys.argv[2],"image":sys.argv[3],"sim_telarray_commit":sys.argv[4],"hessio_rev":sys.argv[5],"host":sys.argv[6]},open(sys.argv[7],"w"),indent=2); open(sys.argv[7],"a").write(chr(10))' \
  "$SIMTOOLS_PROD_IMAGE" "$image_id" "$run_image" "${simtel_commit:-}" "${hessio_rev:-}" "$host" "$output/container.json"
# Same layout as the development container: checkouts read-only, results writable.
models=${VALIDATION_MODELS_DIR:-simulation-models-dev}
mounts=(--volume "$PWD/obdeect-validation:/validation:ro"
        --volume "$PWD/sim_telarray:/reference:ro"
        --volume "$(cd "$models" && pwd):/models:ro")
# CPUS=0 pins the container to one core (podman --cpuset-cpus; 0,1 for two); unset: all cores.
# CASE=cases/performance/air-shower.json runs that one case instead of the suite.
if [[ -n "${CASE:-}" ]]; then
  select=(--case "/validation/$CASE")
else
  select=(--suite "/validation/suites/$suite.json")
fi
podman run --rm "${mounts[@]}" ${extra_env[@]+"${extra_env[@]}"} \
  ${CPUS:+--cpuset-cpus "$CPUS"} \
  --volume "$output:/results:Z" \
  "$run_image" \
  obdeect-validate run \
    "${select[@]}" \
    --environment "/validation/environments/$environment.json" \
    --output /results --non-interactive
