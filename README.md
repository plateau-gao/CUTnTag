# CUTnTag
This is for CUT&Tag data process.
* The `input.csv` contains all the fastq files gathering from NGS. There are five key words serving as column names which cannot be changed: `group`, `cond`, `SampleId`,`read1` and `read2`
* The parameters are set in the `cutntag.config` file. The syntax of nextflow config file is `key=value`. Users can change the 'value' in the config file to change the parameters used in pipeline. Please do not change the 'key', which may cause errors to the pipelines.
* The `cutntag.def` is used for build singularity image for this pipelines.
* For the pipeline it self:
  - Data pre-processing:
    - Quality control of fastq file: `fastqc`(optional)
    - Trim adaptors: `trim-galore` (optional)
    - Build reference genome index: `bowtie2` (optional)
  - Alignment to Hg38 & Spike-in genome
    - Alignment: `bowtie2`
    - Remove duplicates: `picard`. (optional)
    - Assess mapped fragment size
    - Alignment results filtering (optional)
  - File format conversion
    - Sam to Bam file: `samtools`
    - Bam to bed file: `bedtools`
    - Bed to bedgraph file: `bedtools`
    - Bam to bigwig file: `deeptools` (optional, if users want to draw heatmap with `deeptools`)
    - Access reproducibility
  - Peak calling
    - `SEACR`
    - `macs2`
    - Noted: if users do not have control group, the pipeline will use SEACR for peak calling by default, and the numeric threshold is 0.1
  - Visualization (optional)
    - `deeptools`
  - Data summary
    - `R`

> reference: Zheng Y et al (2020). Protocol.io
