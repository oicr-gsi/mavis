# mavis

MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data.

## Overview

## Dependencies

* [mavis 2.2.6](http://mavis.bcgsc.ca/)
* [bcftools 1.9 ](https://samtools.github.io/bcftools/bcftools.html)


## Usage

### Cromwell
```
java -jar cromwell.jar run mavis.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`sampleId`|String|sample identifier, which will be used for final naming of output files
`inputBAMs`|Array[BamData]|Collection of alignment files with indexes and metadata
`svData`|Array[SvData]|Collection of SV calls with metadata
`reference`|String|The genome reference build. for example: hg19, hg38


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filterDellyInput.jobMemory`|Int|12|Memory allocated for this job
`filterDellyInput.timeout`|Int|5|Timeout in hours, needed to override imposed limits
`filterDellyInput.svFileBase`|String|basename(svFile,".vcf.gz")|the basename of the file, needed to form the output file name
`runMavis.outputCONFIG`|String|"mavis_config.cfg"|name of config file for MAVIS
`runMavis.scriptName`|String|"mavis_config.sh"|name for bash script to run mavis configuration, default mavis_config.sh
`runMavis.mavisAligner`|String|"blat"|blat by default, may be customized
`runMavis.mavisScheduler`|String|"SGE"|Our cluster environment, sge, SLURM etc. 
`runMavis.mavisDrawFusionOnly`|String|"False"|flag for MAVIS visualization control
`runMavis.mavisAnnotationMemory`|Int|32000|Memory allocated for annotation step
`runMavis.mavisValidationMemory`|Int|32000|Memory allocated for validation step
`runMavis.mavisTransValidationMemory`|Int|32000|Memory allocated for transvalidation step
`runMavis.mavisMemoryLimit`|Int|32000|Max Memory allocated for MAVIS
`runMavis.mavisQueue`|String|"u20.q"|the mavis job queue
`runMavis.minClusterPerFile`|Int|10|Determines the way parallel calculations are organized 
`runMavis.drawNonSynonymousCdnaOnly`|String|"False"|flag for MAVIS visualization control
`runMavis.mavisUninformativeFilter`|String|"True"|Should be enabled if used is only interested in events inside genes, speeds up calculations
`runMavis.jobMemory`|Int|12|Memory allocated for this job
`runMavis.sleepInterval`|Int|20|A pause after scheduling step, in seconds
`runMavis.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`runMavis.maxBins`|Int|100000|Maximum value for transcriptome_bins and genome_bins parameters, Default is 100000
`runMavis.mavisMaxTime`|Int|timeout * 1800|Timeout for MAVIS tasks, in seconds. 1/2 of the timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`summary`|File|File with copy number variants, native varscan format|vidarr_label: summary
`drawings`|File|Plots generated with MAVIS, collected into a single tar.gz archive|vidarr_label: drawings
`nscvWT`|File?|Whole transcriptome non-synonymous coding variants. The output file is only generated if variants are found|vidarr_label: nscvWT
`nscvWG`|File?|Whole genome non-synonymous coding variants. The output file is only generated if variants are found|vidarr_label: nscvWG


## Commands
This section lists command(s) run by WORKFLOW workflow
 
* Running mavis
 
MAVIS annotates structural variants for WG and WT experiments
 
 
### OPTIONAL : Filter Delly files to keep ONLY the PASS calls
 
```
     bcftools view -i "%FILTER='PASS'" ~{svFile} -Oz -o ~{svFileBase}.pass.vcf.gz
```
 
 
### Setup Mavis : Inline python code
 
```
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
```
 
### Run Mavis
 
```
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
```
 
Compile results.  Drawings and Legends are collected into a single zip archive.  
 
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
