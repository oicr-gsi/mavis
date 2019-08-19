version 1.0

workflow MavisWorkflow {
input {
 	String projectID
	String outputDIR
	String outputCONFIG
        String mavisModule
	File inputBAM
	File inputBAMindex
	File STARFusion
}

call configureMavis { input: outputDIR = outputDIR, outputCONFIG = outputCONFIG, projectID = projectID, STARFusion = STARFusion, inputBAM = inputBAM, inputBAMidx = inputBAMindex, modules = mavisModule }
call setupMavis { input: inputAnchor = configureMavis.output_message, outputDIR = outputDIR, outputCONFIGfile = configureMavis.mavis_config, STARFusion = STARFusion, modules = mavisModule }
call launchMavis { input: batchID = setupMavis.batchID, outputDIR = outputDIR, outputCONFIG = outputCONFIG, modules = mavisModule, projectID = projectID }
}

# ===================================
#            CONFIGURE
# ===================================
task configureMavis {
input {
 	File   inputBAM
 	File   inputBAMidx
	File   STARFusion
        Int?   jobMemory = 10
        File?  referenceGenome = "/.mounts/labs/gsiprojects/gsi/reference/mavis/hg19/hg19_random.fa"
        File?  annotations = "/.mounts/labs/gsiprojects/gsi/reference/mavis/hg19/ensembl69_hg19_annotations_with_ncrna.json"
        File?  masking = "/.mounts/labs/gsiprojects/gsi/reference/mavis/hg19/hg19_masking.tab"
        File?  dvgAnnotations = "/.mounts/labs/gsiprojects/gsi/reference/mavis/hg19/dgv_hg19_variants.tab"
        File?  alignerReference = "/.mounts/labs/gsiprojects/gsi/reference/mavis/hg19/hg19.2bit"
        File?  templateMetadata = "/.mounts/labs/gsiprojects/gsi/reference/mavis/hg19/cytoBand.txt"
	String outputDIR
	String outputCONFIG
	String projectID
        String? modules = "mavis/2.2.6"
}

command <<<
 set -euo pipefail
 export MAVIS_REFERENCE_GENOME=~{referenceGenome}
 export MAVIS_ANNOTATIONS=~{annotations}
 export MAVIS_MASKING=~{masking}
 export MAVIS_DGV_ANNOTATION=~{dvgAnnotations}
 export MAVIS_ALIGNER_REFERENCE=~{alignerReference}
 export MAVIS_TEMPLATE_METADATA=~{templateMetadata}
 mavis config --write "~{outputDIR}~{outputCONFIG}" \
              --assign ~{projectID} starfusion starfusion \
              --library ~{projectID} transcriptome diseased True ~{inputBAM} \
              --convert starfusion ~{STARFusion} starfusion
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  String output_message = read_string(stdout())
  File mavis_config = "${outputDIR}${outputCONFIG}"
}
}


# ====================================
#          SETUP
# ====================================
task setupMavis {
input {
	String outputDIR
	String outputCONFIGfile
	File   STARFusion
        Int?   jobMemory = 10
	String inputAnchor
        String? mavisAligner = "blat"
	String? mavisScheduler = "SGE"
	String? mavisDrawFusionOnly = "False"
	Int? mavisAnnotationMemory = 32000
	Int? mavisValidationMemory = 32000
	Int? mavisTransValidationMemory = 32000
	Int? mavisMemoryLimit = 32000
        Int? minClusterPerFile = 5
	String? drawNonSynonymousCdnaOnly = "False"
        String? mavisUninformativeFilter = "True"        
        String? modules = "mavis/2.2.6"
}

command <<<
  set -euo pipefail
  export MAVIS_ALIGNER='~{mavisAligner}'
  export MAVIS_SCHEDULER=~{mavisScheduler}
  export MAVIS_DRAW_FUSIONS_ONLY=~{mavisDrawFusionOnly}
  export MAVIS_ANNOTATION_MEMORY=~{mavisAnnotationMemory}
  export MAVIS_VALIDATION_MEMORY=~{mavisValidationMemory}
  export MAVIS_TRANS_VALIDATION_MEMORY=~{mavisTransValidationMemory}
  export MAVIS_MEMORY_LIMIT=~{mavisMemoryLimit}
  export DRAW_NON_SYNONYMOUS_CDNA_ONLY=~{drawNonSynonymousCdnaOnly}
  export min_clusters_per_file=~{minClusterPerFile}
  export MAVIS_UNINFORMATIVE_FILTER=~{mavisUninformativeFilter}
  mavis setup ~{outputCONFIGfile} -o ~{outputDIR}
  BATCHID=$(grep MS_batch ~{outputDIR}/build.cfg | grep -v \] | sed s/.*-// | tail -n 1)
  echo $BATCHID
>>>

runtime {
  memory: "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
 String batchID = read_string(stdout())
}

}

# ====================================
#          LAUNCH
# ====================================
task launchMavis {
input {
	String outputDIR
        String projectID
	String outputCONFIG
	String batchID
        String? modules = "mavis/2.2.6"
        Int?   jobMemory = 12
        Int?   sleepInterval = 10
}


command <<<
 echo "launching MAVIS..."
 mavis schedule -o ~{outputDIR} --submit 2> >(tee launch_stderr.log)
 echo "Scheduler launched"
 sleep 10
 LASTJOB=$(cat launch_stderr.log | grep SUBMITTED | tail -n 1 | sed s/.*\(//)
 num='([0-9^]+)'
 if [[ $LASTJOB =~ $num ]]; then
    jobID=$BASH_REMATCH
    echo "Job id: $jobID"
    while qstat | grep $jobID; do
        echo "Mavis Pipeline is still running"
        sleep 5
    done
    if [ -f ~{outputDIR}/summary/MAVIS-$jobID.COMPLETE ]; then
        zip -j ~{batchID}_drawings.zip ~{outputDIR}/~{projectID}\_diseased_transcriptome/annotate/*/drawings/*svg \
                                                ~{outputDIR}/~{projectID}\_diseased_transcriptome/annotate/*/drawings/*json
        exit 0
    fi
    echo "MAVIS job finished but THERE ARE NO RESULTS"
 fi
 exit 1
>>>

runtime {
  memory: "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
 File zipped_drawings = "${batchID}_drawings.zip"
 File mavis_results   = "${outputDIR}/summary/mavis_summary_all_${projectID}.tab"
}

}

