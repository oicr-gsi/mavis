version 1.0

workflow mavis {
input {
 	String donor
        Array[Pair[String, File]] inputBAMs
        Array[Pair[Pair[String, String],File]] svData
        Array[File]   inputBAMidx
}

scatter (i in range(length(inputBAMs))) {
    Int idx = i
    call massageBamData{input: in = inputBAMs[idx]}
}
Array[File] bamFiles = massageBamData.bam
Array[String] libTypes = massageBamData.libType

scatter (j in range(length(svData))) {
    Int jdx = j
    call massageSvData{input: inMetaData = svData[jdx].left, inFile = svData[jdx].right}
}

Array[File] svFiles       = massageSvData.svFile
Array[String] svWorkflows  = massageSvData.svWorkflow
Array[String] svLibDesigns = massageSvData.libDesign

call runMavis { input: donor = donor, svData = svFiles, inputBAMs = bamFiles, inputBAMidx = inputBAMidx, libTypes = libTypes, svWorkflows = svWorkflows, svLibDesigns = svLibDesigns }

meta {
 author: "Peter Ruzanov"
 email: "peter.ruzanov@oicr.on.ca"
 description: "mavis 1.0"
}

output {
  File zippedSummaryTable = select_first(runMavis.zipped_summaries)
  File zippedDrawings = select_first(runMavis.zipped_drawings)
}
}

# ===================================
#   MASSAGE alignment data
# ===================================
task massageBamData{
input{
  Pair[String, File] in
}

parameter_meta {
 in: "paired library design and the corresponding bam file"
}

command <<<
 echo "Processing pair ~{in.left} and  ~{basename(in.right)}"
>>>

output{
 String libType = "~{in.left}"
 File bam = in.right
}
}

# ===================================
#   MASSAGE SV data
# ===================================
task massageSvData{
input{
  Pair[String, String] inMetaData
  File inFile
}

parameter_meta {
 inMetaData: "paired library design and workflow name for SV data file"
 inFile: "SV data file"
}

command <<<
 echo "Processing pair ~{inMetaData.left} and  ~{inMetaData.right}"
 echo "Have file ~{basename(inFile)}"
>>>

output{
 String svWorkflow = "~{inMetaData.left}"
 String libDesign  = "~{inMetaData.right}"
 File svFile = inFile
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
        Array[String] svLibDesigns
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

parameter_meta {
 inputBAMs: "array of input bam files"
 inputBAMidx: "array of input .bai files"
 svData: "array of SV calls"
 libTypes: "List of library types, metadata for inputBAMs"
 svWorkflows: "List of SV callers, metadata for svData"
 svLibDesigns: "List of library designs to accompany the list of SV calls"
 outputCONFIG: "name of config file for MAVIS"
 scriptName: "name for bash script to run mavis configuration, default mavis_config.sh"
 donor: "donor id, i.e. PCSI_0001 Identifies a patient, cell culture grown at certain condition etc."
 referenceGenome: "path to fasta file with genomic assembly"
 annotations: ".json file with annotations for MAVIS"
 masking: "masking data in .tab format"
 dvgAnnotations: "The DGV annotations help to deal with variants found in normal tissue"
 alignerReference: "References in 2bit (compressed) format, used by MAVIS aligner"
 templateMetadata: "Chromosome Band Information, used for visualization"
 mavisAligner: "blat by default, may be customized"
 mavisScheduler: "Our cluster environment, sge, SLURM etc. "
 mavisDrawFusionOnly: "flag for MAVIS visualization control"
 mavisAnnotationMemory: "Memory allocated for annotation step"
 mavisValidationMemory: "Memory allocated for validation step"
 mavisTransValidationMemory: "Memory allocated for transvalidation step"
 mavisMemoryLimit: "Max Memory allocated for MAVIS"
 minClusterPerFile: "Determines the way parallel calculations are organized "
 drawNonSynonymousCdnaOnly: "flag for MAVIS visualization control"
 mavisUninformativeFilter: "Should be enabled if used is only interested in events inside genes, speeds up calculations"
 modules: "modules needed to run MAVIS"
 jobMemory: "Memory allocated for this job"
 sleepInterval: "A pause after scheduling step, in seconds"
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

 b = "~{sep=' ' inputBAMs}"
 bams = b.split()
 l = "~{sep=' ' libTypes}"
 libs = l.split()
 s = "~{sep=' ' svData}"
 svdata = s.split()
 w = "~{sep=' ' svWorkflows}"
 wfs = w.split()
 sl = "~{sep=' ' svLibDesigns}"
 svlibs = sl.split()

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
           assign_arrays[svlibs[s]].append(wfMappings[w])

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

