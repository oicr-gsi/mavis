version 1.0

workflow mavis {
  input {
    String sampleId
    Array[BamData] inputBAMs
    Array[SvData] svData
	String reference
  }

  parameter_meta {
    sampleId: "sample identifier, which will be used for final naming of output files"
    inputBAMs: "Collection of alignment files with indexes and metadata"
    svData: "Collection of SV calls with metadata"
	reference: "The genome reference build. for example: hg19, hg38"
  }
  
  
  String filter_modules = "bcftools/1.9"
  
  Map[String,String] mavis_modules_by_genome = { "hg19": "mavis/2.2.6 mavis-config/1.2 hg19-mavis/2.2.6 hg19/p13", "hg38" : "mavis/2.2.6 mavis-config/1.2 hg38v110-mavis/2.2.6 hg38/p12" }
  String mavis_modules = mavis_modules_by_genome [ reference ]
  
  Map[String,String] resources = { 
  "hg38_annotations": "$HG38V110_MAVIS_ROOT/ensembl_v110_hg38_annotations.json", 
  "hg38_dvgAnnotations": "$HG38_MAVIS_ROOT/dgv_hg38_variants.tab",
  "hg38_cytoband": "$HG38_MAVIS_ROOT/cytoBand.txt",
  "hg38_masking": "$HG38_MAVIS_ROOT/hg38_masking.tab",
  "hg38_referenceGenome": "$HG38_ROOT/hg38_random.fa",
  "hg38_alignerReference": "$HG38_MAVIS_ROOT/hg38.2bit",
  "hg19_annotations": "$HG19_MAVIS_ROOT/ensembl69_hg19_annotations_with_ncrna.json", 
  "hg19_dvgAnnotations": "$HG19_MAVIS_ROOT/dgv_hg19_variants.tab",
  "hg19_cytoband": "$HG19_MAVIS_ROOT/cytoBand.txt",
  "hg19_masking": "$HG19_MAVIS_ROOT/hg19_masking.tab",
  "hg19_referenceGenome": "$HG19_ROOT/hg19_random.fa",
  "hg19_alignerReference": "$HG19_MAVIS_ROOT/hg19.2bit"
 }

  String build_annotations = resources [ "~{reference + '_annotations'}" ]

  String sanitized_sid = sub(sampleId, "_", ".")

  scatter(b in inputBAMs) {
    File bams = b.bam
    File bamIndexes = b.bamIndex
    String bamLibraryDesigns = b.libraryDesign
  }

  scatter(s in svData) {
    if( select_first([s.doFilter,false]) && s.workflowName == "delly"){
      call filterDellyInput {
        input:
          svFile = s.svFile,
          modules = filter_modules
      }
    }
    File svFiles = select_first([filterDellyInput.fsvFile,s.svFile])
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
      svLibDesigns = svLibraryDesigns,
      modules = mavis_modules,
      annotations       = resources [ "~{reference + '_annotations'}" ],
      dvgAnnotations    = resources [ "~{reference + '_dvgAnnotations'}" ],
      templateMetadata  = resources [ "~{reference + '_cytoband'}" ],
      masking           = resources [ "~{reference + '_masking'}" ],
      referenceGenome   = resources [ "~{reference + '_referenceGenome'}" ],
      alignerReference  = resources [ "~{reference + '_alignerReference'}" ],
      arribaConverter   = "$MAVIS_CONFIG_ROOT/bin/parse_arriba.py"
  }

  meta {
   author: "Peter Ruzanov, Lawrence Heisler"
   email: "peter.ruzanov@oicr.on.ca, lawrence.heisler@oicr.on.ca"
   description: "MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data."
   dependencies: [
      {
        name: "mavis/2.2.6",
        url: "http://mavis.bcgsc.ca/"
      },
      {
        name: "bcftools/1.9 ",
        url: "https://samtools.github.io/bcftools/bcftools.html"
      }
    ]
    output_meta: {
    summary: {
        description: "File with copy number variants, native varscan format",
        vidarr_label: "summary"
    },
    drawings: {
        description: "Plots generated with MAVIS, collected into a single tar.gz archive",
        vidarr_label: "drawings"
    },
    nscvWT: {
        description: "Whole transcriptome non-synonymous coding variants. The output file is only generated if variants are found",
        vidarr_label: "nscvWT"
    },
    nscvWG: {
        description: "Whole genome non-synonymous coding variants. The output file is only generated if variants are found",
        vidarr_label: "nscvWG"
    }
}
  }



  output {
    File summary = runMavis.summary
    File drawings = runMavis.drawings
    File? nscvWT   = runMavis.nscvWT
    File? nscvWG   = runMavis.nscvWG
  }
}


# ===================================
# OPTIONAL Filter Input Task
# ===================================
task filterDellyInput {
  input {
    File svFile
    String modules
    Int jobMemory = 12
    Int timeout = 5
    String svFileBase = basename(svFile,".vcf.gz")
  }
  parameter_meta {
    svFile: "the file that needs to be filtered"
    svFileBase: "the basename of the file, needed to form the output file name"
    modules: "modules needed to run filtering"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    bcftools view -i "%FILTER='PASS'" ~{svFile} -Oz -o ~{svFileBase}.pass.vcf.gz
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File fsvFile = "~{svFileBase}.pass.vcf.gz"
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
    String mavisQueue = "u20.q"
    Int minClusterPerFile = 10
    String drawNonSynonymousCdnaOnly = "False"
    String mavisUninformativeFilter = "True"
    String modules
    Int jobMemory = 12
    Int sleepInterval = 20
    Int timeout = 24
    Int maxBins = 100000
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
    mavisQueue: "the mavis job queue"
    minClusterPerFile: "Determines the way parallel calculations are organized "
    drawNonSynonymousCdnaOnly: "flag for MAVIS visualization control"
    mavisUninformativeFilter: "Should be enabled if used is only interested in events inside genes, speeds up calculations"
    modules: "modules needed to run MAVIS"
    jobMemory: "Memory allocated for this job"
    sleepInterval: "A pause after scheduling step, in seconds"
    timeout: "Timeout in hours, needed to override imposed limits"
    maxBins: "Maximum value for transcriptome_bins and genome_bins parameters, Default is 100000"
    mavisMaxTime: "Timeout for MAVIS tasks, in seconds. 1/2 of the timeout"
  }

  command <<<
 
    export MAVIS_REFERENCE_GENOME=~{referenceGenome}
    export MAVIS_ANNOTATIONS=~{annotations}
    export MAVIS_MASKING=~{masking}
    export MAVIS_DGV_ANNOTATION=~{dvgAnnotations}
    export MAVIS_ALIGNER_REFERENCE=~{alignerReference}
    export MAVIS_TEMPLATE_METADATA=~{templateMetadata}
    export MAVIS_TIME_LIMIT=~{mavisMaxTime}
    # we're using system python3
    python3 <<CODE

    libtypes = {'WT': "transcriptome", 'MR': "transcriptome", 'WG': "genome"}
    wfMappings = {'StructuralVariation': 'delly', 'delly': 'delly', 'arriba' : 'arriba', 'starFusion': 'starfusion', 'StarFusion': 'starfusion', 'starfusion': 'starfusion', 'manta': 'manta'}

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

    config_lines = []
    assign_arrays = {}
    for lt in libtypes.keys():
       assign_arrays[lt] = []

    for b in range(len(bams)):
       flag = ('False' if libs[b] == 'WG' else 'True')
       config_lines.append( "--library " + libs[b] + ".~{sid} " + libtypes[libs[b]] + " diseased " + flag + " " + bams[b] + " \\\\" )


    for s in range(len(svdata)):
       for w in wfMappings.keys():
          if w in wfs[s]:
             if w == 'arriba':
                config_lines.append( "--external_conversion arriba \"~{arribaConverter}  " + svdata[s] + "\"" + " \\\\" )
             else:
                config_lines.append( "--convert " + wfMappings[w] + " " + svdata[s] + " " + wfMappings[w] + " \\\\" )
             assign_arrays[svlibs[s]].append(wfMappings[w])

    for b in range(len(bams)):
       if len(assign_arrays[libs[b]]) > 0:
          separator = " "
          tools = separator.join(assign_arrays[libs[b]])
          config_lines.append( "--assign " + libs[b] + ".~{sid} " + tools + " \\\\" )

    f = open("~{scriptName}","w+")
    f.write("#!/bin/bash" + "\n\n")
    f.write('mavis config \\\\\n')
    f.write('\n'.join(config_lines) + '\n')
    if "WT" in libs or "MR" in libs:
       f.write("--transcriptome_bins 500" + ' \\\\\n')
    if "WG" in libs:
       f.write("--genome_bins 500" + ' \\\\\n')
    f.write("--write ~{outputCONFIG}\n")
    f.close()
    CODE
    
    chmod +x ~{scriptName}
    ./~{scriptName} &
    wait
    
    if [ ! -f ~{outputCONFIG} ]; then
      sed -i 's/_bins 500/_bins ~{maxBins}/' ~{scriptName}
      ./~{scriptName}
    fi

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
    export MAVIS_QUEUE=~{mavisQueue}
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
          ### create an empty zip file, which will be updated with drawings and legends.  if there are none, than the empty file is provisioned out
          echo | zip -q > ~{prefix}.mavis_drawings.zip && zip -dq ~{prefix}.mavis_drawings.zip -

          ### find all drawing directories, recursively add the drawings
          for draw_dir in `ls -d *~{sid}\_diseased_*/annotate/*/drawings`
          do
            zip -qjur ~{prefix}.mavis_drawings.zip $draw_dir
          done

          ### there should be a single mavis_summary_all files
          cp summary/mavis_summary_all_*.tab ~{prefix}.mavis_summary.tab

          ### non-synonymous coding variants are separate into WG or WT files; each may or may not be produced
          if [ -e summary/mavis_summary_WG.*_non-synonymous_coding_variants.tab ];then
            cp summary/mavis_summary_WG.*_non-synonymous_coding_variants.tab ~{prefix}.WG_non-synonymous_coding_variants.tab
          fi
          if [ -e summary/mavis_summary_WT.*_non-synonymous_coding_variants.tab ];then
            cp summary/mavis_summary_WT.*_non-synonymous_coding_variants.tab ~{prefix}.WT_non-synonymous_coding_variants.tab
          fi		  
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
    File summary   = "~{prefix}.mavis_summary.tab"
    File? nscvWT   = "~{prefix}.WT_non-synonymous_coding_variants.tab"
    File? nscvWG   = "~{prefix}.WG_non-synonymous_coding_variants.tab"
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
  Boolean? doFilter
}
