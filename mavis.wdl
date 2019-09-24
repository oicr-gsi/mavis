version 1.0

workflow mavis {
input {
 	String projectID
	String? outputCONFIG = "mavis_config.cfg"
        String? mavisModules = "mavis/2.2.6 hg19-mavis/2.2.6 hg19/p13"
	File inputBAM
	File inputBAMindex
	File fusionData
}

call runMavis { input: outputCONFIG = outputCONFIG, projectID = projectID, fusionData = fusionData, inputBAM = inputBAM, inputBAMidx = inputBAMindex, modules = mavisModules }

output {
  File resultTable    = runMavis.mavis_results
  File zippedDrawings = select_first(runMavis.zipped_drawings)
}
}

# ===================================
#  CONFIGURE, SETUP and LAUNCH
# ===================================
task runMavis {
input {
 	File   inputBAM
 	File   inputBAMidx
	File   fusionData
        String? outputCONFIG
        String projectID
        String?  referenceGenome = "$HG19_ROOT/hg19_random.fa"
        String?  annotations = "$HG19_MAVIS_ROOT/ensembl69_hg19_annotations.json"
        String?  masking = "$HG19_MAVIS_ROOT/hg19_masking.tab"
        String?  dvgAnnotations = "$HG19_MAVIS_ROOT/dgv_hg19_variants.tab"
        String?  alignerReference = "$HG19_MAVIS_ROOT/hg19.2bit"
        String?  templateMetadata = "$HG19_MAVIS_ROOT/cytoBand.txt"
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
        String? modules = "mavis/2.2.6 hg19-mavis/2.2.6 hg19/p13"
        Int?   jobMemory = 12
        Int?   sleepInterval = 20
}

command <<<
 unset LD_LIBRARY_PATH
 unset LD_LIBRARY_PATH_modshare
 export MAVIS_REFERENCE_GENOME=~{referenceGenome}
 export MAVIS_ANNOTATIONS=~{annotations}
 export MAVIS_MASKING=~{masking}
 export MAVIS_DGV_ANNOTATION=~{dvgAnnotations}
 export MAVIS_ALIGNER_REFERENCE=~{alignerReference}
 export MAVIS_TEMPLATE_METADATA=~{templateMetadata}
 mavis config --write ~{outputCONFIG} \
              --assign ~{projectID} starfusion starfusion \
              --library ~{projectID} transcriptome diseased True ~{inputBAM} \
              --convert starfusion ~{fusionData} starfusion
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
 mavis setup ~{outputCONFIG} -o .
 BATCHID=$(grep MS_batch build.cfg | grep -v \] | sed s/.*-// | tail -n 1)
 mavis schedule -o . --submit 2> >(tee launch_stderr.log)
 sleep ~{sleepInterval}
 LASTJOB=$(cat launch_stderr.log | grep SUBMITTED | tail -n 1 | sed s/.*\(//)
 num='([0-9^]+)'
 if [[ $LASTJOB =~ $num ]]; then
    jobID=$BASH_REMATCH
    while qstat | grep $jobID; do
        sleep 5
    done
    if [ -f summary/MAVIS-$jobID.COMPLETE ]; then
        zip -qj $BATCHID"_drawings.zip" ~{projectID}\_diseased_transcriptome/annotate/*/drawings/*svg \
                                        ~{projectID}\_diseased_transcriptome/annotate/*/drawings/*json
        exit 0
    fi
    echo "MAVIS job finished but THERE ARE NO RESULTS"
    exit 1
 fi
 echo "Could not retrieve last job id"
 exit 1
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  Array[File?] zipped_drawings = glob('*.zip')
  File mavis_results   = "summary/mavis_summary_all_${projectID}.tab"
}
}

