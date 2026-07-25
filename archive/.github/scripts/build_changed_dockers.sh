#!/usr/bin/env bash
set -euo pipefail

# Builds and pushes every dockerfiles/Dockerfile.* changed in the most recent
# commit, resolving each to its canonical image name via dockerfiles/image_names.env.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MAPPING_FILE="${REPO_ROOT}/dockerfiles/image_names.env"

changed=$(git -C "${REPO_ROOT}" diff --name-only HEAD~1 HEAD -- 'dockerfiles/Dockerfile.*')

if [ -z "${changed}" ]; then
    echo "No Dockerfile changes detected."
    exit 0
fi

while IFS= read -r path; do
    dockerfile_name="$(basename "${path}")"
    suffix="${dockerfile_name#Dockerfile.}"
    image_name=$(grep -E "^${suffix}=" "${MAPPING_FILE}" | cut -d= -f2)

    if [ -z "${image_name}" ]; then
        echo "ERROR: no image name mapping found for ${dockerfile_name} in ${MAPPING_FILE}" >&2
        exit 1
    fi

    echo "---------------------------------------------------"
    echo "Building and pushing ${dockerfile_name} as ${image_name}"
    "${REPO_ROOT}/dockerfiles/build_push.sh" "${image_name}" "${dockerfile_name}"
done <<< "${changed}"
