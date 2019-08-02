#!/bin/bash
set -eu
set -o pipefail

CROMWELL_HOST=""
CROMWELL_WORKFLOW_ID_PATH=""
OUTPUT_DIR=""
OUTPUT_FILE_DEFINITION=""

while (( "$#" )); do
    [[ -z "${2+not_set}" ]] && echo "Missing value for ${1}" && exit 1
    case "${1}" in
		--cromwell-host) CROMWELL_HOST="${2}"; shift 2 ;;
		--cromwell-workflow-id-path) CROMWELL_WORKFLOW_ID_PATH="${2}"; shift 2 ;;
		--output-dir) OUTPUT_DIR="${2}"; shift 2 ;;
		--output-file-definition) OUTPUT_FILE_DEFINITION="${2}"; shift 2 ;;
		--*) echo "Unsupported option ${1}" >&2 ; exit 1 ;;
		*) break ;;
    esac
done

[[ -z "${CROMWELL_HOST}" ]] && echo "--cromwell-host is not set" >&2 && exit 1
[[ -z "${CROMWELL_WORKFLOW_ID_PATH}" ]] && echo "--cromwell-workflow-id-path is not set" >&2 && exit 1
[[ -z "${OUTPUT_DIR}" ]] && echo "--output-dir is not set" >&2 && exit 1
[[ -z "${OUTPUT_FILE_DEFINITION}" ]] && echo "--output-file-definition is not set" >&2 && exit 1

WORKFLOW_ID=$(cat "${CROMWELL_WORKFLOW_ID_PATH}")
if [[ -z "${WORKFLOW_ID}" ]]; then
  echo "Unable to read workflow id file." >&2
  exit 1
fi

STATUS_TEXT=$(curl -s -X GET "${CROMWELL_HOST}/api/workflows/v1/${WORKFLOW_ID}/status")
if [[ $(jq --raw-output --exit-status '.status' <<< "${STATUS_TEXT}") != "Succeeded" ]]; then
  echo "Status response from cromwell:" >&2
  echo "${STATUS_TEXT}" >&2
fi

OUTPUTS_TEXT=$(curl -s -X GET "${CROMWELL_HOST}/api/workflows/v1/${WORKFLOW_ID}/outputs")
if ! jq --raw-output --exit-status '.outputs' <<< "${OUTPUTS_TEXT}" >/dev/null ; then
  echo "Output response from cromwell:" >&2
  echo "${OUTPUTS_TEXT}" >&2
else
  OUTPUT=$(jq --raw-output --exit-status '.outputs' <<< "${OUTPUTS_TEXT}")
fi

echo "${OUTPUT}"
echo "Implementation not complete" >&2
exit 1
