version 1.0

workflow mavis {
  input {
    String donor
    Array[BamData] inputBAMs
    Array[SvData] svData
    Map[String, String] libTypeMap = {"WT": "transcriptome", "MR": "transcriptome", "WG": "genome"}
  }

  String sanitized_donor = sub(donor, "_", ".")
  String configFileName = "mavis_config.cfg"
  String scriptFileName = "mavis_config.sh"

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

  call generateConfigScript {
    input:
    donor = sanitized_donor,
    inputBAMs = bams,
    inputBAMidx = bamIndexes,
    libTypes = bamLibraryDesigns,
    libTypeMap = libTypeMap,
    svData = svFiles,
    svWorkflows = workflowNames,
    svLibDesigns = svLibraryDesigns,
    configFileName = configFileName,
    scriptFileName = scriptFileName
  }

  call config {
    input:
    configScript = generateConfigScript.configScript,
    configName = configFileName
  }

  call setup {
    input:
    configFile = config.configFile
  }

  scatter(ct in setup.clusterTabs) {
    call validate {
      input: inFiles=ct, script=setup.submitValidate
    }
    call annotate {
      input: inFiles=validate.outFile, script=setup.submitAnnotate
    }
  }

  call submitMultiple as pairing {
    input: inFiles=annotate.outFile, script=setup.submitPairing, outDirName="pairing"
  }

  call submitMultiple as summary {
    input: inFiles=pairing.outFiles, script=setup.submitSummary, outDirName="summary"
  }

  call zipResults{
    #input: inFiles=flatten([summary.outFiles, flatten(annotate.drawings)]), batchID=setup.batchID, donor=sanitized_donor
	input: inFiles=flatten([flatten(annotate.drawings)]), batchID=setup.batchID, donor=sanitized_donor
  }

  output {
    File results = zipResults.zipArchive
    Array[File] summaries = summary.outFiles
  }

  meta {
   author: "Peter Ruzanov, Iain Bancarz, Lawrence Heisler"
   email: "peter.ruzanov@oicr.on.ca, ibancarz@oicr.on.ca, lawrence.heisler@oicr.on.ca"
   description: "MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data. The workflow runs each MAVIS action as a WDL task; this replaces the MAVIS default method, of submitting actions directly to a computing cluster."
   dependencies: [
      {
        name: "mavis/2.2.6",
        url: "http://mavis.bcgsc.ca/"
      }
    ]
    output_meta: {
      summaryFiles: "File with copy number variants, native varscan format",
      zippedDrawings: "File of plots generated with MAVIS"
    }
  }

  parameter_meta {
    donor: "Donor id"
    inputBAMs: "Collection of alignment files with indexes and metadata"
    svData: "Collection of SV calls with metadata"
    libTypeMap: "Mapping from library type strings to Mavis library types"
  }

}

task generateConfigScript {

  meta {
    description: "Python code to generate the Mavis configuration script"
    output_meta: {
      configScript: "Bash script for Mavis configuration"
    }
  }

  parameter_meta {
    inputBAMs: "Array of input BAM files"
    inputBAMidx: "Array of input BAM index files"
    svData: "Array of somatic variant data files"
    libTypes: "Array of library type strings"
    libTypeMap: "Mapping from library type strings to Mavis library types"
    svWorkflows: "Array of somatic variant workflow strings"
    svLibDesigns: "Array of somatic variant library design strings"
    configFileName: "Name of file output from the config bash script"
    scriptFileName: "Name of the config bash script"
    arribaConverter: "Path to arriba conversion script"
    donor: "String representing the donor"
    jobMemory: "Memory for the task, in gigabytes"
    modules: "Environment modules for the task"
    timeout: "Timeout for the task, in hours"
  }

  input {
    Array[File] inputBAMs
    Array[File] inputBAMidx
    Array[File] svData
    Array[String] libTypes
    Map[String, String] libTypeMap
    Array[String] svWorkflows
    Array[String] svLibDesigns
    String configFileName = "mavis_config.cfg"
    String scriptFileName = "mavis_config.sh"
    String arribaConverter
    String donor
    String modules
    Int jobMemory = 12
    Int timeout = 24
  }
    
  command <<<
    python <<CODE
    import json

    libtypes = json.loads(open("~{write_json(libTypeMap)}").read())
    wfMappings = {'StructuralVariation': 'delly', 'delly': 'delly', 'arriba' : 'arriba', 'StarFusion': 'starfusion', 'starFusion': 'starfusion', 'manta': 'manta'}

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
             if w == 'arriba':
                 convert_lines.append( "--external_conversion arriba \"~{arribaConverter}  " + svdata[s] + "\"" + " \\\\" )
             else:
                 convert_lines.append( "--convert " + wfMappings[w] + " " + svdata[s] + " " + wfMappings[w] + " \\\\" )
             assign_arrays[svlibs[s]].append(wfMappings[w])

    for b in range(len(bams)):
       if len(assign_arrays[libs[b]]) > 0:
           separator = " "
           tools = separator.join(assign_arrays[libs[b]])
           assign_lines.append( "--assign " + libs[b] + ".~{donor} " + tools + " \\\\" )

    f = open("~{scriptFileName}","w+")
    f.write("#!/bin/bash" + "\n\n")
    f.write('mavis config \\\\\n')
    f.write('\n'.join(library_lines) + '\n')
    f.write('\n'.join(convert_lines) + '\n')
    f.write('\n'.join(assign_lines) + '\n')
    f.write("--write ~{configFileName}\n")
    f.close()
    CODE
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File configScript = scriptFileName
  }

}

task config {

  meta {
    description: "Run the Mavis config script, writing a config file for setup."
    output_meta: {
      configFile: "Output file written by the config script"
    }
  }

  parameter_meta {
    configScript: "Configuration Bash script for Mavis"
    configName: "Name of output configuration file"
    referenceGenome: "Path to reference FASTA file"
    annotations: "Path to annotations JSON file for MAVIS"
    masking: "Path to reference masking file in .tab format"
    dgvAnnotations: "Path to reference Database of Genome Variants (DGV) annotations"
    alignerReference: "Path to reference in 2bit (compressed) format, for the Mavis aligner"
    templateMetadata: "Chromosome Band Information, used for visualization"
    jobMemory: "Memory for the task, in gigabytes"
    modules: "Environment modules for the task"
    timeout: "Timeout for the task, in hours"
  }

  input {
    File configScript
    String configName
    String referenceGenome
    String annotations
    String masking
    String dgvAnnotations
    String alignerReference
    String templateMetadata
    String modules
    Int jobMemory = 12
    Int timeout = 24
  }

  # needs environment variables -- will run successfully without them, but leave blank params which cause failure in next task

  command <<<
    set -euo pipefail
    export MAVIS_REFERENCE_GENOME=~{referenceGenome}
    export MAVIS_ANNOTATIONS=~{annotations}
    export MAVIS_MASKING=~{masking}
    export MAVIS_DGV_ANNOTATION=~{dgvAnnotations}
    export MAVIS_ALIGNER_REFERENCE=~{alignerReference}
    export MAVIS_TEMPLATE_METADATA=~{templateMetadata}
    chmod +x ~{configScript}
    ~{configScript}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File configFile = configName
  }
}

task setup {

  meta {
    description: "Run Mavis setup; clusters inputs for parallelization, and generates submit.sh scripts for the Mavis pipeline. Also parses the batch ID."
    output_meta: {
      clusterTabs: "Array of TSV files, with input data split for parallelization",
      submitValidate: "submit.sh script for validation",
      submitAnnotate: "submit.sh script for annotation",
      submitPairing: "submit.sh script for pairing",
      submitSummary: "submit.sh script for summary",
      batchID: "Batch ID string"
    }
  }

  parameter_meta {
    configFile: "Mavis config file"
    mavisAligner: "String identifying the aligner"
    mavisDrawFusionOnly: "Mavis parameter: Must be 'True' or 'False'"
    minClusterPerFile: "Minimum number of clusters per file"
    drawNonSynonymousCdnaOnly: "Mavis parameter: Must be 'True' or 'False'"
    mavisUninformativeFilter: "Mavis parameter: Must be 'True' or 'False'"
    mavisMaxFiles: "Maximum number of chunks for parallelization"
    jobMemory: "Memory for the task, in gigabytes"
    modules: "Environment modules for the task"
    timeout: "Timeout for the task, in hours"
  }

  input {
    File configFile
    String mavisAligner = "blat"
    String mavisDrawFusionOnly = "False"
    Int minClusterPerFile = 5
    String drawNonSynonymousCdnaOnly = "False"
    String mavisUninformativeFilter = "True"
    Int mavisMaxFiles=200
    String modules
    Int jobMemory = 32
    Int timeout = 24
  }

  String batchFileName = "batch.txt"

  command <<<
    set -euo pipefail
    unset LD_LIBRARY_PATH
    unset LD_LIBRARY_PATH_modshare
    export MAVIS_ALIGNER='~{mavisAligner}'
    export MAVIS_DRAW_FUSIONS_ONLY=~{mavisDrawFusionOnly}
    export MAVIS_MAX_FILES=~{mavisMaxFiles}
    export DRAW_NON_SYNONYMOUS_CDNA_ONLY=~{drawNonSynonymousCdnaOnly}
    export min_clusters_per_file=~{minClusterPerFile}
    export MAVIS_UNINFORMATIVE_FILTER=~{mavisUninformativeFilter}
    mavis setup ~{configFile} -o .
    # Try to get batch ID from the "build.cfg" file generated by setup
    grep MS_batch build.cfg | grep -v \] | sed s/.*-// | tail -n 1 | sed s/\n// > ~{batchFileName} || \
      echo "unknown_batch" > ~{batchFileName}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Array[File] clusterTabs = glob("*/cluster/batch*.tab")
    File submitValidate = glob("*/validate/submit.sh")[0]
    File submitAnnotate = glob("*/annotate/submit.sh")[0]
    File submitPairing = "pairing/submit.sh"
    File submitSummary = "summary/submit.sh"
    String batchID = read_string("~{batchFileName}")
  }
}


task validate {

  # edit and run the validate submit.sh script generated by Mavis

  meta {
    description: "Edit and run the validation submit.sh script generated by Mavis. Editing enables script to be run as a WDL task, instead of submission to a cluster by `mavis schedule`."
    output_meta: {
      outFile: "TSV file with data which passed validation"
    }
  }

  parameter_meta {
    inFiles: "Space-separated list of input files"
    script: "submit.sh script to edit and run"
    jobMemory: "Memory for the task, in gigabytes"
    modules: "Environment modules for the task"
    timeout: "Timeout for the task, in hours"
  }

  input {
    String inFiles # space-separated list
    File script
    String modules
    Int jobMemory = 32
    Int timeout = 24
  }

  command <<<
    set -euo pipefail
    sed -i "s|^[[:space:]]*cd .*||g" ~{script} # do not change the working directory
    sed -i "s|^[[:space:]]*echo .*||g" ~{script} # do not echo start/finish times
    sed -i "s|--inputs .*|--inputs ~{inFiles} \\\\|g" ~{script}
    sed -i "s|--output .*|--output .|g" ~{script}
    chmod +x ~{script}

    ~{script}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File outFile = "validation-passed.tab"
  }
}

task annotate {

  # edit and run the annotate submit.sh script generated by Mavis
  # similar to validate, but also return contents of the drawings directory (if any)

  meta {
    description: "Edit and run the annotation submit.sh script generated by Mavis. Editing enables script to be run as a WDL task, instead of submission to a cluster by `mavis schedule`. Used for the Mavis pairing and summary steps."
    output_meta: {
      outFile: "TSV file with annotated data",
      drawings: "Files with .json and .svg plot output"
    }
  }

  parameter_meta {
    inFiles: "Space-separated list of input files"
    script: "submit.sh script to edit and run"
    jobMemory: "Memory for the task, in gigabytes"
    modules: "Environment modules for the task"
    timeout: "Timeout for the task, in hours"
  }

  input {
    String inFiles # space-separated list
    File script
    String modules
    Int jobMemory = 32
    Int timeout = 24
  }

  command <<<
    set -euo pipefail
    sed -i "s|^[[:space:]]*cd .*||g" ~{script} # do not change the working directory
    sed -i "s|^[[:space:]]*echo .*||g" ~{script} # do not echo start/finish times
    sed -i "s|--output .*|--output . \\\\|g" ~{script}
    sed -i "s|--inputs .*|--inputs ~{inFiles}|g" ~{script}
    chmod +x ~{script}
    ~{script}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File outFile = "annotations.tab"
    Array[File] drawings = glob("drawings/{*json,*svg}")
  }
}

task submitMultiple {

  # edit and run a submit.sh script generated by Mavis, with multiple inputs

  meta {
    description: "Edit and run a submit.sh script generated by Mavis with multiple inputs and outputs. Editing enables script to be run as a WDL task, instead of submission to a cluster by `mavis schedule`"
    output_meta: {
      outFiles: "TSV files with output data"
    }
  }

  parameter_meta {
    inFiles: "Space-separated list of input files"
    script: "submit.sh script to edit and run"
    outDirName: "Name of output directory; will be created as a subdirectory of the task working directory"
    jobMemory: "Memory for the task, in gigabytes"
    modules: "Environment modules for the task"
    timeout: "Timeout for the task, in hours"
  }

  input {
    Array[File] inFiles
    File script
    String outDirName
    String modules
    Int jobMemory = 32
    Int timeout = 24
  }

  # order of --inputs and --output is inconsistent; always follow with a \

  command <<<
    set -euo pipefail
    mkdir ~{outDirName}
    sed -i "s|^[[:space:]]*cd .*||g" ~{script} # do not change the working directory
    sed -i "s|^[[:space:]]*echo .*||g" ~{script} # do not echo start/finish times
    sed -i "s|--inputs .*|--inputs ~{sep=' ' inFiles} \\\\|g" ~{script}
    sed -i "s|--output .*|--output ~{outDirName} \\\\|g" ~{script}
    chmod +x ~{script}
    ~{script}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Array[File] outFiles = glob("~{outDirName}/*\.tab")
  }
}

task zipResults {

  meta {
    description: "Gather required result files into a .zip archive"
    output_meta: {
      zipArchive: "ZIP archive file"
    }
  }

  parameter_meta {
    inFiles: "Array of input files"
    batchID: "Batch ID string"
    donor: "Donor ID string"
    jobMemory: "Memory for the task, in gigabytes"
    modules: "Environment modules for the task"
    timeout: "Timeout for the task, in hours"
  }

  input {
    Array[File] inFiles
    String batchID
    String donor
    String modules
    Int jobMemory = 12
    Int timeout = 24
  }

  # create a directory for the zip archive; allows unzip without exploding multiple files into the working directory

  #String outPrefix = "~{batchID}.~{donor}_mavis-output"
  String outPrefix = "~{donor}_drawings"

  command <<<
    set -euo pipefail
    mkdir ~{outPrefix}
    #cp -t ~{outPrefix} ~{sep=' ' inFiles}
    for file in ~{sep=' ' inFiles}
    do
      cp $file ~{outPrefix}
    done
    zip -qr ~{outPrefix}.zip ~{outPrefix}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File zipArchive = "~{outPrefix}.zip"
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
