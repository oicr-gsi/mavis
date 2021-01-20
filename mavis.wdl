version 1.0

workflow mavis {
  input {
    String donor
    Array[BamData] inputBAMs
    Array[SvData] svData
    String configFileName = "mavis_config.cfg"
    String scriptFileName = "mavis_config.sh"
  }

  String sanitized_donor = sub(donor, "_", ".")

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

  call generateConfigFile {
    input:
    donor = sanitized_donor,
    inputBAMs = bams,
    inputBAMidx = bamIndexes,
    libTypes = bamLibraryDesigns,
    svData = svFiles,
    svWorkflows = workflowNames,
    svLibDesigns = svLibraryDesigns,
    configFileName = configFileName,
    scriptFileName = scriptFileName
  }

  call config {
    input:
    configScript = generateConfigFile.configScript,
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

  call zipResults {
    input: drawings=flatten(annotate.drawings), summaries=summary.outFiles, batchID=setup.batchID, donor=sanitized_donor
  }

  meta {
   author: "Peter Ruzanov, Iain Bancarz"
   email: "peter.ruzanov@oicr.on.ca, ibancarz@oicr.on.ca"
   description: "MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data."
   dependencies: [
      {
        name: "mavis/2.2.6",
        url: "http://mavis.bcgsc.ca/"
      }
    ]
    output_meta: {
      zippedSummary: "File with copy number variants, native varscan format",
      zippedDrawings: "Plots generated with MAVIS"
    }
  }

  parameter_meta {
    donor: "Donor id"
    inputBAMs: "Collection of alignment files with indexes and metadata"
    svData: "Collection of SV calls with metadata"
  }

  output {
    Array[File?] zippedDrawings = zipResults.zippedDrawings
    Array[File?] zippedSummary = zipResults.zippedSummary
  }
}

task generateConfigFile {

  input {
    Array[File] inputBAMs
    Array[File] inputBAMidx
    Array[File] svData
    Array[String] libTypes
    Array[String] svWorkflows
    Array[String] svLibDesigns
    String configFileName = "mavis_config.cfg"
    String scriptFileName = "mavis_config.sh"
    String arribaConverter
    String donor
    String modules
    Int jobMemory = 12
    Int sleepInterval = 20
    Int timeout = 24
    Int mavisMaxTime = timeout * 1800
  }
    
  command <<<
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

  # needs environment variables -- will run successfully without them, but leave blank params which cause failure in next task

  input {
    File configScript
    String configName
    String referenceGenome
    String annotations
    String masking
    String dvgAnnotations
    String alignerReference
    String templateMetadata
    String modules
    Int jobMemory = 12
    Int sleepInterval = 20
    Int timeout = 24
  }

  command <<<
    set -euo pipefail
    export MAVIS_REFERENCE_GENOME=~{referenceGenome}
    export MAVIS_ANNOTATIONS=~{annotations}
    export MAVIS_MASKING=~{masking}
    export MAVIS_DGV_ANNOTATION=~{dvgAnnotations}
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

  input {
    File configFile
    String mavisAligner = "blat"
    String mavisDrawFusionOnly = "False"
    Int minClusterPerFile = 5
    String drawNonSynonymousCdnaOnly = "False"
    String mavisUninformativeFilter = "True"
    Int mavisMaxFiles=200
    String modules
    Int jobMemory = 12
    Int sleepInterval = 20
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
    # MS_batch may be empty; do not error on a non-zero return code from grep
    grep MS_batch ~{configFile} | grep -v \] | sed s/.*-// | tail -n 1 | sed s/\n// > ~{batchFileName} || true
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

  input {
    String inFiles # space-separated list
    File script
    String modules
    Int jobMemory = 12
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

  output {
    File outFile = "validation-passed.tab"
  }
}

task annotate {

  # edit and run the annotate submit.sh script generated by Mavis
  # similar to validate, but also return contents of the drawings directory (if any)

  input {
    String inFiles # space-separated list
    File script
    String modules
    Int jobMemory = 12
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

  output {
    File outFile = "annotations.tab"
    Array[File] drawings = glob("drawings/{*json,*svg}")
  }
}

task submitMultiple {

  # edit and run a submit.sh script generated by Mavis, with multiple inputs

  input {
    Array[File] inFiles
    File script
    String outDirName
    String modules
    Int jobMemory = 12
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

  output {
    Array[File] outFiles = glob("~{outDirName}/*\.tab")
  }
}

task zipResults {

  input {
    Array[File] drawings
    Array[File] summaries
    String batchID
    String donor
    String modules
    Int jobMemory = 12
    Int timeout = 24
  }

  # if batchID is empty, replace with a default value
  # do not set -euo pipefail; causes glob failure if drawings.zip is absent

  command <<<
    if [ -z ~{batchID} ]; then
      BATCH="unknown-batch"
    else
      BATCH=~{batchID}
    fi
    zip -qj $BATCH.~{donor}_drawings.zip ~{sep=' ' drawings}
    zip -qj $BATCH.~{donor}_summary.zip ~{sep=' ' summaries}
  >>>

  output {
    Array[File?] zippedDrawings = glob('*drawings.zip')
    Array[File?] zippedSummary = glob('*summary.zip')
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
