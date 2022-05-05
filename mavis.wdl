version 1.0

workflow mavis {
  input {
    String sampleId
    Array[BamData] inputBAMs
    Array[SvData] svData
  }

  String sanitized_sid = sub(sampleId, "_", ".")

  scatter(b in inputBAMs) {
    File bams = b.bam
    File bamIndexes = b.bamIndex
    String bamLibraryDesigns = b.libraryDesign
  }

  scatter(s in svData) {
    File svFiles = s.svFile
    String workflowNames = s.workflowName
    String svLibraryDesigns = s.libraryDesign
  }

  call runMavis {
    input:
      sid = sanitized_sid,
      prefix = sampleId,
      inputBAMs = bams,
      inputBAMidx = bamIndexes,
      libTypes = bamLibraryDesigns,
      svData = svFiles,
      svWorkflows = workflowNames,
      svLibDesigns = svLibraryDesigns
  }

  meta {
   author: "Peter Ruzanov, Lawrence Heisler"
   email: "peter.ruzanov@oicr.on.ca, lawrence.heisler@oicr.on.ca"
   description: "MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data."
   dependencies: [
      {
        name: "mavis/2.2.6",
        url: "http://mavis.bcgsc.ca/"
      }
    ]
    output_meta: {
      summary: "File with copy number variants, native varscan format",
      drawings: "Plots generated with MAVIS, collected into a single tar.gz archve"
    }
  }

  parameter_meta {
    sampleId: "sample identifier, which will be used for final naming of output files"
    inputBAMs: "Collection of alignment files with indexes and metadata"
    svData: "Collection of SV calls with metadata"
  }

  output {
    File summary = runMavis.summary
    File drawings = runMavis.drawings
  }
}

# ===================================
#  CONFIGURE, SETUP and LAUNCH
# ===================================
task runMavis {
  input {
    Array[File] inputBAMs
    Array[File] inputBAMidx
    Array[File] svData
    Array[String] libTypes
    Array[String] svWorkflows
    Array[String] svLibDesigns
    String outputCONFIG = "mavis_config.cfg"
    String scriptName = "mavis_config.sh"
    String arribaConverter
    String sid
    String prefix
    String referenceGenome
    String annotations
    String masking
    String dvgAnnotations
    String alignerReference
    String templateMetadata
    String mavisAligner = "blat"
    String mavisScheduler = "SGE"
    String mavisDrawFusionOnly = "False"
    Int mavisAnnotationMemory = 32000
    Int mavisValidationMemory = 32000
    Int mavisTransValidationMemory = 32000
    Int mavisMemoryLimit = 32000
    Int minClusterPerFile = 10
    String drawNonSynonymousCdnaOnly = "False"
    String mavisUninformativeFilter = "False"
    String modules
    Int jobMemory = 12
    Int sleepInterval = 20
    Int timeout = 24
    Int mavisMaxTime = timeout * 1800
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
    arribaConverter: "path to arriba conversion script"
    sid: "sample ID, this is provided to maivs and cannot include reseerved characters [;,_\\s] "
    prefix: "the final sid which will be used for naming of the output files"
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
    timeout: "Timeout in hours, needed to override imposed limits"
    mavisMaxTime: "Timeout for MAVIS tasks, in seconds. 1/2 of the timeout"
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
    export MAVIS_TIME_LIMIT=~{mavisMaxTime}
    python <<CODE

    libtypes = {'WT': "transcriptome", 'MR': "transcriptome", 'WG': "genome"}
    wfMappings = {'StructuralVariation': 'delly', 'delly': 'delly', 'arriba' : 'arriba', 'StarFusion': 'starfusion', 'manta': 'manta'}

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
     library_lines.append( "--library " + libs[b] + ".~{sid} " + libtypes[libs[b]] + " diseased " + flag + " " + bams[b] + " \\\\" )


    for s in range(len(svdata)):
     for w in wfMappings.keys():
         if w in wfs[s]:
             if w == 'arriba':
                 convert_lines.append( "--external_conversion arriba \"~{arribaConverter}  " + svdata[s] + "\"" + " \\\\" )
             else:
                 convert_lines.append( "--convert " + wfMappings[w] + " " + svdata[s] + " " + wfMappings[w] + " \\\\" )
             assign_arrays[svlibs[s]].append(wfMappings[w])

    for b in range(len(bams)):
       if len(assign_arrays[libs[b]]) > 0:
           separator = " "
           tools = separator.join(assign_arrays[libs[b]])
           assign_lines.append( "--assign " + libs[b] + ".~{sid} " + tools + " \\\\" )

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
          zip -qj ~{prefix}".mavis_drawings.zip" *~{sid}\_diseased_*/annotate/*/drawings/*svg \
                                                *~{sid}\_diseased_*/annotate/*/drawings/*json
          cp summary/mavis_summary_all_*~{sid}.tab ~{prefix}.mavis.summary.tab
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
    timeout: "~{timeout}"
  }

  output {
    File drawings  = "~{prefix}.mavis_drawings.zip"
    File summary   = "~{prefix}.mavis.summary.tab"
  }
}

struct BamData {
  File bam
  File bamIndex
  String libraryDesign
}

struct SvData {
  File svFile
  String workflowName
  String libraryDesign
}
