#!/bin/bash
set -eu
set -o pipefail

JAVA=""
SEQWARE_JAR=""
WORKFLOW_RUN_SWID=""
PROCESSING_SWID=""
CROMWELL_HOST=""
CROMWELL_WORKFLOW_ID_PATH=""
OUTPUT_DIR=""

while (( "$#" )); do
    [[ -z "${2+not_set}" ]] && echo "Missing value for ${1}" && exit 1
    case "${1}" in
		--java-path) JAVA="${2}"; shift 2 ;;
		--seqware-jar-path) SEQWARE_JAR="${2}"; shift 2 ;;
		--workflow-run-swid) WORKFLOW_RUN_SWID="${2}"; shift 2 ;;
		--processing-swid) PROCESSING_SWID="${2}"; shift 2 ;;
		--cromwell-host) CROMWELL_HOST="${2}"; shift 2 ;;
		--cromwell-workflow-id-path) CROMWELL_WORKFLOW_ID_PATH="${2}"; shift 2 ;;
		--output-dir) OUTPUT_DIR="${2}"; shift 2 ;;
		--*) echo "Unsupported option ${1}" >&2 ; exit 1 ;;
		*) break ;;
    esac
done

[[ -z "${JAVA}" ]] && echo "--java-path is not set" >&2 && exit 1
[[ -z "${SEQWARE_JAR}" ]] && echo "--seqware-jar-path is not set" >&2 && exit 1
[[ -z "${WORKFLOW_RUN_SWID}" ]] && echo "--workflow-run-swid is not set" >&2 && exit 1
[[ -z "${PROCESSING_SWID}" ]] && echo "--processing-swid is not set" >&2 && exit 1
[[ -z "${CROMWELL_HOST}" ]] && echo "--cromwell-host is not set" >&2 && exit 1
[[ -z "${CROMWELL_WORKFLOW_ID_PATH}" ]] && echo "--cromwell-workflow-id-path is not set" >&2 && exit 1
[[ -z "${OUTPUT_DIR}" ]] && echo "--output-dir is not set" >&2 && exit 1

CROMWELL_WORKFLOW_ID=$(cat "${CROMWELL_WORKFLOW_ID_PATH}")
if [[ -z "${CROMWELL_WORKFLOW_ID}" ]]; then
  echo "Unable to read cromwell workflow id file." >&2
  exit 1
fi

STATUS_TEXT=$(curl -s -X GET "${CROMWELL_HOST}/api/workflows/v1/${CROMWELL_WORKFLOW_ID}/status")
if [[ $(jq --raw-output --exit-status '.status' <<< "${STATUS_TEXT}") != "Succeeded" ]]; then
  echo "Status response from cromwell:" >&2
  echo "${STATUS_TEXT}" >&2
fi

OUTPUTS_TEXT=$(curl -s -X GET "${CROMWELL_HOST}/api/workflows/v1/${CROMWELL_WORKFLOW_ID}/outputs")
if ! jq --raw-output --exit-status '.outputs' <<< "${OUTPUTS_TEXT}" >/dev/null ; then
  echo "Outputs response from cromwell:" >&2
  echo "${OUTPUTS_TEXT}" >&2
  exit 1
else
  # store outputs in an array: output file key, output file path
  OUTPUTS=( $(jq --raw-output --exit-status '.outputs | to_entries[] | "\(.key),\(.value)"' <<< "${OUTPUTS_TEXT}") )
  echo "Detected ${#OUTPUTS[@]} output files"
fi

for i in ${!OUTPUTS[@]}; do
  IFS=, read -ra REC <<< "${OUTPUTS[$i]}"

  # figure out the file metatype
  case "$(basename ${REC[1]})" in
    *.bam)
      METATYPE="application/bam"
    ;;
    *.bai)
      METATYPE="application/bam-index"
    ;;
    *.g.vcf.gz)
      METATYPE="application/g-vcf-gz"
    ;;
    *.json)
      METATYPE="text/json"
    ;;
    *.pdf)
      METATYPE="application/pdf"
    ;;
    *.tar.gz|*.tgz)
      METATYPE="application/tar-gzip"
    ;;
    *.tbi)
      METATYPE="application/tbi"
    ;;
    *.vcf.gz)
      METATYPE="application/vcf-gz"
    ;;
    *.zip)
      METATYPE="application/zip-report-bundle"
    ;;
    *.fastq.gz)
      METATYPE="chemical/seq-na-fastq-gzip"
    ;;
    *.fastq)
      METATYPE="chemical/seq-na-fastq"
    ;;
    *.png)
      METATYPE="image/png"
    ;;
    *.bed|*.BedGraph)
      METATYPE="text/bed"
    ;;
    *.fpkm_tracking)
      METATYPE="text/fpkm-tracking"
    ;;
    *.gtf)
      METATYPE="text/gtf"
    ;;
    *.html)
      METATYPE="text/html"
    ;;
    *.vcf)
      METATYPE="text/vcf"
    ;;
    *.txt.gz|*.gz)
      METATYPE="application/txt-gz"
    ;;
    *.out|*.log|*.txt)
      METATYPE="txt/plain"
    ;;
    *)
      METATYPE="unknown/unknown"
    ;;
  esac

  OUTPUT_PATH="${OUTPUT_DIR}$(basename ${REC[1]})"
  echo "Provisioning out ${REC[0]} (metatype = $METATYPE) to ${OUTPUT_PATH}"

  PROCESSING_FILE="s${WORKFLOW_RUN_SWID}_${PROCESSING_SWID}_provision_out_${i}_processing_accession"
  "${JAVA}" \
  -XX:+UseSerialGC -Xmx500M \
  -classpath "${SEQWARE_JAR}" \
  net.sourceforge.seqware.pipeline.runner.Runner \
  --metadata \
  --metadata-parent-accession "${PROCESSING_SWID}" \
  --metadata-workflow-run-ancestor-accession "${WORKFLOW_RUN_SWID}" \
  --metadata-processing-accession-file "${PROCESSING_FILE}" \
  --metadata-processing-accession-file-lock "${PROCESSING_FILE}.lock" \
  --module net.sourceforge.seqware.pipeline.modules.utilities.ProvisionFiles \
  -- \
  --input-file-metadata "pfo::${METATYPE}::${REC[1]}" \
  --output-file "${OUTPUT_PATH}" \
  --force-copy

done
