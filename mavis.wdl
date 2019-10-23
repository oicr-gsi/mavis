version 1.0

workflow mavis {
input {
 	String donor
        Array[File]   inputBAMs
        Array[File]   inputBAMidx
        Array[File]   svData
        Array[String] libTypes
        Array[String] svWorkflows
}

call runMavis { input: donor = donor, svData = svData, inputBAMs = inputBAMs, inputBAMidx = inputBAMidx, libTypes = libTypes, svWorkflows = svWorkflows }

output {
  File zippedSummaryTable = select_first(runMavis.zipped_summaries)
  File zippedDrawings = select_first(runMavis.zipped_drawings)
}
}

# ===================================
#  CONFIGURE, SETUP and LAUNCH
# ===================================
task runMavis {
input {
 	Array[File]   inputBAMs
 	Array[File]   inputBAMidx
	Array[File]   svData
        Array[String] libTypes
        Array[String] svWorkflows
        String? outputCONFIG = "mavis_config.cfg"
        String? scriptName = "mavis_config.sh"
        String donor
        String?  referenceGenome = "$HG19_ROOT/hg19_random.fa"
        String?  annotations = "$HG19_MAVIS_ROOT/ensembl69_hg19_annotations_with_ncrna.json"
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
 python <<CODE

 libtypes = {'WT': "transcriptome", 'MR': "transcriptome", 'WG': "genome"}
 wfMappings = {'StructuralVariation': 'delly', 'Delly': 'delly', 'StarFusion': 'starfusion', 'Manta': 'manta'}
 svMappings = {'delly': ["WG"], 'starfusion': ["WT","MR"], 'manta': ['WG']}

 b = "~{sep=' ' inputBAMs}"
 bams = b.split()
 l = "~{sep=' ' libTypes}"
 libs = l.split()
 s = "~{sep=' ' svData}"
 svdata = s.split()
 w = "~{sep=' ' svWorkflows}"
 wfs = w.split()

 library_lines = []
 convert_lines = []
 assign_lines = []
 assign_arrays = {}
 for lt in libtypes.keys():
   assign_arrays[lt] = []

 for b in range(len(bams)):
   flag = ('False' if libs[b] == 'WG' else 'True')
   library_lines.append( "--library " + libs[b] + ".~{donor} " + libtypes[libs[b]] + " diseased " + flag + " " + bams[b] + " \\\\" )


 for s in range(len(svdata)):
   for w in wfMappings.keys():
       if w in wfs[s]:
           convert_lines.append( "--convert " + wfMappings[w] + " " + svdata[s] + " " + wfMappings[w] + " \\\\" )
           for library_type in svMappings[wfMappings[w]]:
              assign_arrays[library_type].append(wfMappings[w])

 for b in range(len(bams)):
     if len(assign_arrays[libs[b]]) > 0:
         separator = " "
         tools = separator.join(assign_arrays[libs[b]])
         assign_lines.append( "--assign " + libs[b] + ".~{donor} " + tools + " \\\\" )

 f = open("~{scriptName}","w+")
 f.write("#!/bin/bash" + "\n\n")
 f.write('mavis config \\\\\n')
 f.write('\n'.join(library_lines) + '\n')
 f.write('\n'.join(convert_lines) + '\n')
 f.write('\n'.join(assign_lines) + '\n')
 f.write("--write ~{outputCONFIG}\n")
 f.close()
 CODE
 chmod +x ~{scriptName}
 ./~{scriptName}
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
        zip -qj $BATCHID"_drawings.zip" *~{donor}\_diseased_*/annotate/*/drawings/*svg \
                                        *~{donor}\_diseased_*/annotate/*/drawings/*json
        zip -qj $BATCHID"_summary.zip" summary/mavis_summary_all_*~{donor}.tab
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
  Array[File?] zipped_drawings  = glob('*drawings.zip')
  Array[File?] zipped_summaries = glob('*summary.zip')
}
}

