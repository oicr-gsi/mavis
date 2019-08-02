#!/bin/bash
set -eu
set -o pipefail

JAVA=""
CROMWELL_JAR=""
CROMWELL_HOST=""
WORKFLOW=""
INPUTS=""
CROMWELL_WORKFLOW_ID_PATH=""
OPTIONS=""
DEPS_ZIP=""

while (( "$#" )); do
    [[ -z "${2+not_set}" ]] && echo "Missing value for ${1}" && exit 1
    case "${1}" in
		--java-path) JAVA="${2}"; shift 2 ;;
		--cromwell-jar-path) CROMWELL_JAR="${2}"; shift 2 ;;
		--cromwell-host) CROMWELL_HOST="${2}"; shift 2 ;;
		--workflow) WORKFLOW="${2}"; shift 2 ;;
		--inputs) INPUTS="${2}"; shift 2 ;;
		--cromwell-workflow-id-path) CROMWELL_WORKFLOW_ID_PATH="${2}"; shift 2 ;;
		--options) OPTIONS="${2}"; shift 2 ;;
		--deps-zip) DEPS_ZIP="${2}"; shift 2 ;;
		--*) echo "Unsupported option ${1}" >&2 ; exit 1 ;;
		*) break ;;
    esac
done

[[ -z "${JAVA}" ]] && echo "--java-path is not set" >&2 && exit 1
[[ -z "${CROMWELL_JAR}" ]] && echo "--cromwell-jar-path is not set" >&2 && exit 1
[[ -z "${CROMWELL_HOST}" ]] && echo "--cromwell-host is not set" >&2 && exit 1
[[ -z "${WORKFLOW}" ]] && echo "--workflow is not set" >&2 && exit 1
[[ -z "${INPUTS}" ]] && echo "--inputs is not set" >&2 && exit 1
[[ -z "${CROMWELL_WORKFLOW_ID_PATH}" ]] && echo "--cromwell-workflow-id-path is not set" >&2 && exit 1

SUBMIT_CMD="${JAVA} -XX:+UseSerialGC -Xmx1g -jar ${CROMWELL_JAR} submit ${WORKFLOW} --inputs ${INPUTS} --host ${CROMWELL_HOST}"

if [[ -z ${OPTIONS} ]]; then
  LAUNCH_CMD+=" --options ${OPTIONS}"
fi

if [[ -z ${DEPS_ZIP} ]]; then
  LAUNCH_CMD+=" --imports ${DEPS_ZIP}"
fi

SUBMIT_TEXT=$(eval ${SUBMIT_CMD})
echo "${SUBMIT_TEXT}"

WORKFLOW_ID=$(sed -n 's#^\(.*Workflow \)\(.*\)\( submitted to .*\)$#\2#p' <<< "${SUBMIT_TEXT}")
if [[ -z "${WORKFLOW_ID}" ]]; then
  echo "Unable to parse workflow id." >&2
  exit 1
else
  echo "Workflow id = ${WORKFLOW_ID}"
  echo "${WORKFLOW_ID}" > "${CROMWELL_WORKFLOW_ID_PATH}"
fi

STATUS="Pending"
while [[ "${STATUS}" == "Pending" || "${STATUS}" == "Submitted" || "${STATUS}" == "Running" ]]; do
  echo "Workflow id = ${WORKFLOW_ID} status = ${STATUS}"
  sleep 30
  STATUS=$(curl -s -X GET "${CROMWELL_HOST}/api/workflows/v1/${WORKFLOW_ID}/status" | jq --raw-output --exit-status '.status')
done

echo "Workflow id = ${WORKFLOW_ID} status = ${STATUS}"
if [[ "${STATUS}" != "Succeeded" ]]; then
  exit 1
else
  exit 0
fi
