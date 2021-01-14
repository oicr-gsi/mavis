version 1.0

workflow write_config_draft {
  input {
    String donor
    Array[BamData] inputBAMs
    Array[SvData] svData
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

  # TODO break up runMavis into multiple sub-tasks
  # TODO are annotation/validation memory limits needed if running jobs via WDL?

  call generateMavisConfigFile {
    input:
      donor = sanitized_donor,
      inputBAMs = bams,
      inputBAMidx = bamIndexes,
      libTypes = bamLibraryDesigns,
      svData = svFiles,
      svWorkflows = workflowNames,
      svLibDesigns = svLibraryDesigns
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
      zippedSummaryTable: "File with copy number variants, native varscan format",
      zippedDrawings: "Plots generated with MAVIS"
    }
  }

  parameter_meta {
    donor: "Donor id"
    inputBAMs: "Collection of alignment files with indexes and metadata"
    svData: "Collection of SV calls with metadata"
  }

  output {
    File configScript = generateMavisConfigFile.configScript
  }
}

task generateMavisConfigFile {

  input {
    Array[File] inputBAMs
    Array[File] inputBAMidx
    Array[File] svData
    Array[String] libTypes
    Array[String] svWorkflows
    Array[String] svLibDesigns
    String outputConfig = "mavis_config.cfg"
    String scriptName = "mavis_config.sh"
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

    f = open("~{scriptName}","w+")
    f.write("#!/bin/bash" + "\n\n")
    f.write('mavis config \\\\\n')
    f.write('\n'.join(library_lines) + '\n')
    f.write('\n'.join(convert_lines) + '\n')
    f.write('\n'.join(assign_lines) + '\n')
    f.write("--write ~{outputConfig}\n")
    f.close()
    CODE
  >>>

   runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  
}
  output {
    File configScript = scriptName
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