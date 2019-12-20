# mavis

MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data.

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

## Optional Assembly-specific Parameters:

hg19-specific data, for other assemblies these should be changed:

Paramter|Value
---|---
dvgAnnotations | String? (optional, default $HG19_MAVIS_ROOT/dgv_hg19_variants.tab)
templateMetadata | String? (optional, default $HG19_MAVIS_ROOT/cytoBand.txt)
annotations | String? (optional, default $HG19_MAVIS_ROOT/ensembl69_hg19_annotations_with_ncrna.json)
masking | String? (optional, default $HG19_MAVIS_ROOT/hg19_masking.tab)
referenceGenome | String? (optional, default $HG19_ROOT/hg19_random.fa)
alignerReference | String? (optional, default $HG19_MAVIS_ROOT/hg19.2bit)


## Other Parameters with default values:

Paramter|Value
---|---
mavisValidationMemory | Int? (optional, default = 32000)
mavisTransValidationMemory | Int? (optional, default = 32000)
mavisMemoryLimit | Int? (optional, default = 32000)
minClusterPerFile | Int? (optional, default = 5)
mavisAnnotationMemory | Int? (optional, default = 32000)
mavisUninformativeFilter | String? (optional, default = True)
drawNonSynonymousCdnaOnly | String? (optional, default = False)
sleepInterval | Int? (optional, default = 20)
modules | String? (optional, default = "mavis/2.2.6 hg19-mavis/2.2.6 hg19/p13")
outputCONFIG | String? (optional, default = "mavis_config.cfg")
mavisAligner | String? (optional, default = "blat")
mavisScheduler | String? (optional, default = "SGE")
scriptName | String? (optional, default = "mavis_config.sh")
mavisDrawFusionOnly | String? (optional, default = "False")
jobMemory | Int? (optional, default = 12)

## Required Inputs:

Paramter|Value
---|---
inputBAMs | Array[Pair[String, File]] type (WG, WT, MR) and path to input bam file
inputBAMidx | Array[File] indexes of input bam files
svData | Array[Pair[Pair[String, String], File]] tool name, type (WG, WT, MR) and path to file with AV calls
donor | String, id for the sample. Olive will use donor

## Outputs

```
  zippedDrawings     - all drawings which MAVIS produced
  zippedSummaryTable - all summary tables which MAVIS produced

```

