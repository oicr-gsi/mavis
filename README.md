# mavis

MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data. The workflow runs each MAVIS action as a WDL task; this replaces the MAVIS default method, of submitting actions directly to a computing cluster.

## Overview

## Dependencies

* [mavis 2.2.6](http://mavis.bcgsc.ca/)

## Mavis Pipeline

### Input Formats

The workflow expects inputs for SV calls in one of a few supported formats. Currently they are manta, delly and starfusion (the latter one is for WT). The respective olive will collect whatever data are available and run Mavis.

### Setup

The Mavis pipeline is configured and set up as normal. For example:

```
 mavis config
   --library MY.LIBRARY genome diseased False input.bam
   --convert delly delly.merged.vcf.gz delly
   --assign MY.LIBRARY delly
   --write mavis_config.cfg

 mavis setup mavis_config.cfg -o OUTPUT_DIR
```

See [Mavis documentation](https://mavis.readthedocs.io/en/latest/) for further details.

### Running

The workflow modifies Mavis jobs to run using WDL.

After running `mavis setup`, the conventional Mavis pipeline would continue with `mavis submit`. Using `mavis submit` is not desirable for this workflow, because it submits jobs directly to a computing cluster, which are not tracked by OICR workflow systems.

Therefore, the workflow edits the `submit.sh` scripts generated by `mavis setup` and runs them as WDL tasks.


## Usage

### Cromwell
```
java -jar cromwell.jar run mavis.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`donor`|String|Donor id
`inputBAMs`|Array[BamData]|Collection of alignment files with indexes and metadata
`svData`|Array[SvData]|Collection of SV calls with metadata
`generateConfigScript.arribaConverter`|String|Path to arriba conversion script
`generateConfigScript.modules`|String|Environment modules for the task
`config.referenceGenome`|String|Path to reference FASTA file
`config.annotations`|String|Path to annotations JSON file for MAVIS
`config.masking`|String|Path to reference masking file in .tab format
`config.dgvAnnotations`|String|Path to reference Database of Genome Variants (DGV) annotations
`config.alignerReference`|String|Path to reference in 2bit (compressed) format, for the Mavis aligner
`config.templateMetadata`|String|Chromosome Band Information, used for visualization
`config.modules`|String|Environment modules for the task
`setup.modules`|String|Environment modules for the task
`validate.modules`|String|Environment modules for the task
`annotate.modules`|String|Environment modules for the task
`pairing.modules`|String|Environment modules for the task
`summary.modules`|String|Environment modules for the task
`zipResults.modules`|String|Environment modules for the task

#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`generateConfigScript.jobMemory`|Int|12|Memory for the task, in gigabytes
`generateConfigScript.timeout`|Int|24|Timeout for the task, in hours
`config.jobMemory`|Int|12|Memory for the task, in gigabytes
`config.timeout`|Int|24|Timeout for the task, in hours
`setup.mavisAligner`|String|"blat"|String identifying the aligner
`setup.mavisDrawFusionOnly`|String|"False"|Mavis parameter: Must be 'True' or 'False'
`setup.minClusterPerFile`|Int|5|Minimum number of clusters per file
`setup.drawNonSynonymousCdnaOnly`|String|"False"|Mavis parameter: Must be 'True' or 'False'
`setup.mavisUninformativeFilter`|String|"True"|Mavis parameter: Must be 'True' or 'False'
`setup.mavisMaxFiles`|Int|200|Maximum number of chunks for parallelization
`setup.jobMemory`|Int|32|Memory for the task, in gigabytes
`setup.timeout`|Int|24|Timeout for the task, in hours
`validate.jobMemory`|Int|32|Memory for the task, in gigabytes
`validate.timeout`|Int|24|Timeout for the task, in hours
`annotate.jobMemory`|Int|32|Memory for the task, in gigabytes
`annotate.timeout`|Int|24|Timeout for the task, in hours
`pairing.jobMemory`|Int|32|Memory for the task, in gigabytes
`pairing.timeout`|Int|24|Timeout for the task, in hours
`summary.jobMemory`|Int|32|Memory for the task, in gigabytes
`summary.timeout`|Int|24|Timeout for the task, in hours
`zipResults.jobMemory`|Int|12|Memory for the task, in gigabytes
`zipResults.timeout`|Int|24|Timeout for the task, in hours

### Outputs

Output | Type | Description
---|---|---
`zippedSummary`|File|File with copy number variants, native varscan format
`zippedDrawings`|File?|File of plots generated with MAVIS


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
