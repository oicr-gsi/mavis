version 1.0

workflow MavisWorkflow {
input {
 	String projectID
	String outputDIR
	String outputCONFIG
        String pythonModule
	File inputBAM
	File inputBAMindex
	File STARFusion
}

call configureMavis { input: outputDIR = outputDIR, outputCONFIG = outputCONFIG, projectID = projectID, STARFusion = STARFusion, inputBAM = inputBAM, inputBAMidx = inputBAMindex, modules = pythonModule }
call setupMavis { input: inputAnchor = configureMavis.output_message, outputDIR = outputDIR, outputCONFIGfile = configureMavis.mavis_config, STARFusion = STARFusion, modules = pythonModule }
call launchMavis { input: inputAnchor = setupMavis.output_message, outputDIR = outputDIR, outputCONFIG = outputCONFIG, modules = pythonModule }
call getJobID { input: inputStrings = launchMavis.output_message }
call finishMavis { input: lastJobId = getJobID.job_id, outputDIR = outputDIR, projectID = projectID }
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
        File?  referenceGenome = "/.mounts/labs/TGL/gsi/tools/mavis/reference_inputs/hg19_random.fa"
        File?  annotations = "/.mounts/labs/TGL/gsi/tools/mavis/reference_inputs/ensembl69_hg19_annotations_with_ncrna.json"
        File?  masking = "/.mounts/labs/TGL/gsi/tools/mavis/reference_inputs/hg19_masking.tab"
        File?  dvgAnnotations = "/.mounts/labs/TGL/gsi/tools/mavis/reference_inputs/dgv_hg19_variants.tab"
        File?  alignerReference = "/.mounts/labs/TGL/gsi/tools/mavis/reference_inputs/hg19.2bit"
        File?  templateMetadata = "/.mounts/labs/TGL/gsi/tools/mavis/reference_inputs/cytoBand.txt"
	String outputDIR
	String outputCONFIG
	String projectID
        String? modules = "python-gsi/3.6.4"
}

command <<<
 set -euo pipefail
 module load ~{modules}
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
        String? modules = "python-gsi/3.6.4"
}

command <<<
  set -euo pipefail
  module load ~{modules}
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
>>>

runtime {
  memory: "~{jobMemory} GB"
}

output {
 String output_message = read_string(stdout())
}

}

# ====================================
#          LAUNCH
# ====================================
task launchMavis {
input {
	String outputDIR
	String outputCONFIG
	String inputAnchor
        String? modules = "python-gsi/3.6.4"
        Int?   jobMemory = 12
}


command <<<
 set -euo pipefail
 echo "launching MAVIS..."
 module load ~{modules}
 mavis schedule -o ~{outputDIR} --submit
>>>

runtime {
  memory: "~{jobMemory} GB"
  continueOnReturnCode: [0, 1]
}

output {
 Array[String] output_message = read_lines(stderr())
}

}

# ===========================================
#              GET ID
# ===========================================
task getJobID {
input {
  Array[String] inputStrings
}

command <<<
     cat ~{write_lines(inputStrings)} | grep SUBMITTED | tail -n 1 | sed s/.*\(// | sed s/\).*//
>>>

output {
     Int job_id = read_int(stdout())
}
}

# ===========================================
#              FINISH
# ===========================================
task finishMavis {
  input {
	Int    lastJobId
	Int?   sleepInterval = 10
        Int?   jobMemory = 6
	String outputDIR
        String projectID
  }

command <<<
 while qstat -j ~{lastJobId}; do
  sleep ~{sleepInterval}
 done
 if [ -f ~{outputDIR}summary/MAVIS-~{lastJobId}.COMPLETE ]
    then
     zip drawings.zip ~{outputDIR}~{projectID}_diseased_transcriptome/annotate/*/drawings/*svg \
                  ~{outputDIR}~{projectID}_diseased_transcriptome/annotate/*/drawings/*json
     exit 0;
 fi
 exit 1;
>>>

runtime {
  memory:  "~{jobMemory} GB"
}

output {
  File zipped_drawings = "drawings.zip"
  File mavis_results   = "${outputDIR}/summary/mavis_summary_all_${projectID}.tab"
}
}

