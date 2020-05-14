# mavis

MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data.

## Overview

## Dependencies

* [mavis 2.2.6](http://mavis.bcgsc.ca/)

## Usage

## Cromwell

``` 
 java -jar cromwell.jar run mavis.wdl --inputs inputs.json 

```

## Running Pipeline

```
 mavis config
   --library MY.LIBRARY genome diseased False inpput.bam
   --convert delly delly.merged.vcf.gz delly
   --assign MY.LIBRARY delly
   --write mavis_config.cfg

 mavis setup mavis_config.cfg -o OUTPUT_DIR
 mavis schedule -o OUTPUT_DIR
```

The workflow will expect that inputs for SV calls will fall into a few supported formats. Currently they are manta, delly and starfusion (the latter one is for WT). The respective olive will collect whatever data are available and 

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`donor`|String|Donor id
`inputBAMs`|Array[BamData]|Collection of alignment files with indexes and metadata
`svData`|Array[SvData]|Collection of SV calls with metadata
`runMavis.referenceGenome`|String|path to fasta file with genomic assembly
`runMavis.annotations`|String|.json file with annotations for MAVIS
`runMavis.masking`|String|masking data in .tab format
`runMavis.dvgAnnotations`|String|The DGV annotations help to deal with variants found in normal tissue
`runMavis.alignerReference`|String|References in 2bit (compressed) format, used by MAVIS aligner
`runMavis.templateMetadata`|String|Chromosome Band Information, used for visualization
`runMavis.modules`|String|modules needed to run MAVIS


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`runMavis.outputCONFIG`|String|"mavis_config.cfg"|name of config file for MAVIS
`runMavis.scriptName`|String|"mavis_config.sh"|name for bash script to run mavis configuration, default mavis_config.sh
`runMavis.mavisAligner`|String|"blat"|blat by default, may be customized
`runMavis.mavisScheduler`|String|"SGE"|Our cluster environment, sge, SLURM etc.
`runMavis.mavisDrawFusionOnly`|String|"False"|flag for MAVIS visualization control
`runMavis.mavisAnnotationMemory`|Int|32000|Memory allocated for annotation step
`runMavis.mavisValidationMemory`|Int|32000|Memory allocated for validation step
`runMavis.mavisTransValidationMemory`|Int|32000|Memory allocated for transvalidation step
`runMavis.mavisMemoryLimit`|Int|32000|Max Memory allocated for MAVIS
`runMavis.minClusterPerFile`|Int|5|Determines the way parallel calculations are organized
`runMavis.drawNonSynonymousCdnaOnly`|String|"False"|flag for MAVIS visualization control
`runMavis.mavisUninformativeFilter`|String|"True"|Should be enabled if used is only interested in events inside genes, speeds up calculations
`runMavis.jobMemory`|Int|12|Memory allocated for this job
`runMavis.sleepInterval`|Int|20|A pause after scheduling step, in seconds
`runMavis.timeout`|Int|24|Timeout in hours, needed to override imposed limits


### Outputs

Output | Type | Description
---|---|---
`zippedSummaryTable`|File?|File with copy number variants, native varscan format
`zippedDrawings`|File?|Plots generated with MAVIS


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```
* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
