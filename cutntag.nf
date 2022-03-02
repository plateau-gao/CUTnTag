/*
=======================================================================================
    Workflow recording
=======================================================================================
*/
process workflow_record {
    publishDir "${params.result}", mode: 'copy'

    output:
    path("workflow.txt")

    script:
    """
    echo "Whether use control group: ${params.use_ctrl}" >> workflow.txt
    echo "do_fastqc: ${params.workflow.fastqc}" >> workflow.txt
    echo "do_trim: ${params.workflow.trim}" >> workflow.txt
    echo "do_bowtie2_build_spikeIn: ${params.workflow.bowtie2_build_ecoli}" >> workflow.txt
    echo "do_bowtie2_build_hg38: ${params.workflow.bowtie2_build_hg38}" >> workflow.txt
    echo "do_picard_rmDup: ${params.workflow.picard_rmDup}" >> workflow.txt
    echo "do_bam2bigwig: ${params.workflow.bam2bigwig}" >> workflow.txt
    echo "do_map_quality_filter: ${params.workflow.map_quality_filter}" >> workflow.txt
    echo "do_peak_calling_by_seacr: ${params.workflow.seacr}" >> workflow.txt
    echo "do_peak_calling_by_macs2: ${params.workflow.macs2}" >> workflow.txt
    echo "draw heatmap on transcription units: ${params.workflow.heatmap_on_trans}" >> workflow.txt
    echo "draw heatmap on peaks: ${params.workflow.heatmap_on_peak}" >> workflow.txt
    """
}

/*
=======================================================================================
    Data Pre-Processing
=======================================================================================
*/

Channel
    .fromPath(params.input)
    .splitCsv(header:true)
    .map {row -> tuple(row.group, row.cond, row.SampleId, file("${params.data}/${row.SampleId}_${row.read1}"), file("${params.data}/${row.SampleId}_${row.read2}"))}
    .into {fastqc_ch; align_ch}

Channel.fromPath(params.ref_hg38).set{ref_hg38_ch}

Channel.fromPath(params.ref_ecoli).set{ref_ecoli_ch}

if (params.workflow.fastqc) {
    process fastqc {
        publishDir "${params.result}/fastqc", mode:"copy"

        input:
        tuple val(group), val(cond), val(id), path(read1), path(read2) from fastqc_ch

        output:
        path "*"

        script:
        """
        fastqc -f fastq ${read1}
        fastqc -f fastq ${read2}
        """
    }
}

if (params.workflow.trim) {
    process trim {
        publishDir "${params.result}/trimed_fastqa", pattern: "*.fq.gz", mode:"copy"
        publishDir "${params.result}/trimed_fastqa/trimmed_repoort", pattern: "*.txt", mode:"copy"

        input:
        tuple val(group), val(cond), val(id), path(read1), path(read2) from align_ch

        output:
        tuple val(group), val(cond), val(id), path("*1.fq.gz"), path("*2.fq.gz") into (hg38_ch, ecoli_ch)
        path ("*.txt")

        script:
        """
        trim_galore ${params.args.trim_galore} --paired ${read1} ${read2}
        """
    }
} else {align_ch.into {hg38_ch; ecoli_ch}}

if (params.workflow.bowtie2_build_hg38) {
    process build_index_hg38{
        publishDir "${params.result}/index_hg38", mode:"copy"

        input:
        path ref_hg38 from ref_hg38_ch

        output:
        path "hg38.*" into index_hg38_ch

        script:
        """
        bowtie2-build --threads ${params.threads} ${ref_hg38} hg38
        """
    }
} else {ref_hg38_ch.set{index_hg38_ch}}

if (params.workflow.bowtie2_build_ecoli) {
    process build_index_ecoli{
        publishDir "${params.result}/index_ecoli", mode:"copy"

        input:
        path ref_ecoli from ref_ecoli_ch

        output:
        path "ecoli.*" into index_ecoli_ch

        script:
        """
        bowtie2-build --threads ${params.threads} ${ref_ecoli} ecoli
        """    
    }
} else {ref_ecoli_ch.set{index_ecoli_ch}}

/*
=======================================================================================
    Alignment to Hg38 and spike-in genome
=======================================================================================
*/

hg38_ch
    .combine(index_hg38_ch.collect())
    .map {group, cond, id, read1,read2, bt1, bt2, bt3, bt4, bt5, bt6 ->
        [group, cond, id, read1, read2, [bt1,bt2,bt3,bt4,bt5,bt6]]
    }
    .set {align_hg38_ch}

process align_hg38 {
    publishDir "${params.result}/alignment/hg38_sam", mode:"copy"

    input:
    tuple val(group), val(cond), val(id), path(read1), path(read2), path(index) from align_hg38_ch

    output:
    tuple val(group), val(cond), val(id), path("${group}_${id}_bowtie2.sam") into (sam_hg38_picard_ch, sam_hg38_nopicard_ch)
    path "${group}_${id}_bowtie2_hg38.txt" into bowtie2_report_hg38_ch

    script:
    index_base = index[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
    """
    bowtie2 ${params.args.bowtie2_align_target} -p ${params.threads} \
    -x $index_base -1 ${read1} -2 ${read2} \
    -S ${group}_${id}_bowtie2.sam &>${group}_${id}_bowtie2_hg38.txt
    """
}

ecoli_ch
    .combine(index_ecoli_ch.collect())
    .map {group, cond, id, read1,read2, bt1, bt2, bt3, bt4, bt5, bt6 ->
        [group, cond, id, read1, read2, [bt1,bt2,bt3,bt4,bt5,bt6]]
    }
    .set {align_ecoli_ch}

process align_ecoli {
    publishDir "${params.result}/alignment/spikeIn_sam", mode:"copy"

    input:
    tuple val(group), val(cond), val(id), path(read1), path(read2), path(index) from align_ecoli_ch

    output:
    tuple val(group), val(cond), val(id), path("${group}_${id}_bowtie2_spikeIn.sam")
    path "${group}_${id}_bowtie2_spikeIn.txt" into bowtie2_report_spikeIn_ch
    tuple val(group), val(cond), val(id), path ("${group}_${id}_spikeIn.seqDepth") into (seqdepth_4_bedgraph_ch, seqdepth_4_R_ch)

    script:
    index_base = index[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
    """
    # align to ecoli
    bowtie2 ${params.args.bowtie2_align_spikein} -p ${params.threads} \
    -x $index_base -1 ${read1} -2 ${read2} \
    -S ${group}_${id}_bowtie2_spikeIn.sam &> ${group}_${id}_bowtie2_spikeIn.txt

    # count seqDepth
    seqDepthDouble=`samtools view ${params.args.samtool_seqDepth_view} ${group}_${id}_bowtie2_spikeIn.sam | wc -l`
    seqDepth=\$((seqDepthDouble/2))
    echo \$seqDepth > ${group}_${id}_spikeIn.seqDepth
    """
}

/*
=======================================================================================
    Remove Duplicates
=======================================================================================
*/

process picard_rmDup {
    publishDir "${params.result}/alignment/hg38_sam", pattern:"*_sorted.sam", mode:"copy"
    publishDir "${params.result}/alignment/removeDuplicate", pattern:"*_*Dup.sam", mode:"copy"
    publishDir "${params.result}/alignment/removeDuplicate/picard_summary", pattern:"*_*Dup.txt", mode:"copy"

    input:
    tuple val(group), val(cond), val(id), path(ori_sam) from sam_hg38_picard_ch

    output:
    tuple val(group), val(cond), val(id), path("${group}_${id}_sorted.sam")
    tuple val(group), val(cond), val(id), path("${group}_${id}_markDup.sam")
    path ("${group}_${id}_markDup.txt")
    tuple val(group), val(cond), val(id), path("${group}_${id}_rmDup.sam") into rmdup_sam_ch
    path "${group}_${id}_rmDup.txt" into rmdup_report_ch

    script:
    """
    picard SortSam -I ${ori_sam} -O ${group}_${id}_sorted.sam ${params.args.picard_sort}
    picard MarkDuplicates -I ${group}_${id}_sorted.sam -O ${group}_${id}_markDup.sam --METRICS_FILE ${group}_${id}_markDup.txt
    picard MarkDuplicates -I ${group}_${id}_sorted.sam -O ${group}_${id}_rmDup.sam --REMOVE_DUPLICATES true -METRICS_FILE ${group}_${id}_rmDup.txt
    """
}

if (params.workflow.picard_rmDup) {
    rmdup_sam_ch.into{sam4filter_ch; sam4assess_frag}
} else {sam_hg38_nopicard_ch.into{sam4filter_ch; sam4assess_frag}}

/*
=======================================================================================
    Assess Mapped Fragment Size
=======================================================================================
*/

process assess_frag_size_distribution {
    publishDir "${params.result}/alignment/sam/fragmentLen", mode: 'copy'

    input:
    tuple val(group), val(cond), val(id), path(sam) from sam4assess_frag

    output:
    path("${group}_${id}_fragmentLen.txt") into fragmentLen_ch

    script:
    """
    samtools view ${params.args.samtool_assess_frag} ${sam} | ${params.args.assess_frag} > ${group}_${id}_fragmentLen.txt
    """
}

/*
=======================================================================================
    Alignment results filtering
=======================================================================================
*/

if (params.args.filter_miniQualityScore > 0) {params.workflow.map_quality_filter = true}

if (params.workflow.map_quality_filter) {
    process quality_filter {
        publishDir "${params.result}/alignment/hg38_filterScore_${params.agrs.filter_miniQualityScore}_sam", mode:"copy"

        input:
        tuple val(group), val(cond), val(id), path(sam) from sam4filter_ch

        output:
        tuple val(group), val(cond), val(id), path("${group}_${id}_*.sam") into sam2bam_ch

        script:
        """
        samtools view -q ${params.agrs.filter_miniQualityScore} ${sam} >${group}_${id}_qualityScore${params.agrs.filter_miniQualityScore}.sam
        """
    }
} else {sam4filter_ch.set{sam2bam_ch}}

/*
=======================================================================================
    File format conversion
=======================================================================================
*/

process sam2bam {
    publishDir "${params.result}/alignment/bam", mode: "copy"

    input:
    tuple val(group), val(cond), val(id), path(sam) from sam2bam_ch

    output:
    tuple val(group), val(cond), val(id), path("${group}_${id}_mapped.bam") into (bam2bed_ch, bam2bigwig_ch, bam4macs2_ch, bam4seacr_report_ch, bam4macs2_report_ch)
    tuple val(group), val(cond), val(id), path("${group}_${id}_sorted.bam") into sort_bam_ch

    script:
    """
    samtools view ${params.args.samtool_sam2bam_view} ${sam} > ${group}_${id}_mapped.bam
    samtools sort ${group}_${id}_mapped.bam -o ${group}_${id}_sorted.bam
    """
}


process bam2bed {
    publishDir "${params.result}/alignment/bed", mode: "copy"

    input:
    tuple val(group), val(cond), val(id), path(mapped_bam) from bam2bed_ch

    output:
    tuple val(group), val(cond), val(id), path("${group}_${id}.bed") into bed_ch
    tuple val(group), val(cond), val(id), path("${group}_${id}_clean.bed") into clean_bed_ch
    tuple val(group), val(cond), val(id), path("${group}_${id}_fragments.bed") into (fragBed_assess_ch, fragBed_convert_ch)
    
    script:
    """
    bedtools bamtobed -i ${mapped_bam} ${params.args.bedtool_bam2bed} > ${group}_${id}.bed
    awk ${params.args.clean_bed} ${group}_${id}.bed > ${group}_${id}_clean.bed
    cut ${params.args.frag_bed_cut} ${group}_${id}_clean.bed | sort ${params.args.frag_bed_sort} > ${group}_${id}_fragments.bed
    """
}

seqdepth_4_bedgraph_ch
    .join(fragBed_convert_ch, by:[0,1,2])
    .set{bed2bedgraph_ch}
    
process bed2bedgraph {
    echo true
    publishDir "${params.result}/alignment/bedgraph", mode: "copy"

    input:
    tuple val(group), val(cond), val(id), path(seq_Depth), path(fragbed) from bed2bedgraph_ch

    output:
    tuple val(group), val(cond), val(id), path ("${group}_${id}_*.bedgraph") into bedgraph_ch

    script:
    """
    seqDepth=`cat ${seq_Depth}`
    if [[ "\$seqDepth" -gt "1" ]]; then
        scale_factor=`echo "${params.args.constant_C} / \$seqDepth"| bc -l`
        bedtools genomecov -bg -scale \$scale_factor -i ${fragbed} -g ${params.genome_size} > ${group}_${id}_fragments.normalized.bedgraph
    fi
    """
}

if (params.workflow.heatmap_on_peak || params.workflow.heatmap_on_trans) {params.workflow.bam2bigwig = true}

if (params.workflow.bam2bigwig) {

    sort_bam_ch
        .filter {group, cond, id, path -> cond == 'exp'}
        .set{exp_bam_ch}

    process bam2bigwig {
        publishDir "${params.result}/alignment/bigwig", mode:"copy"

        input:
        tuple val(group), val(cond), val(id), path(sort_bam) from exp_bam_ch

        output:
        tuple val(group), val(id), path("${group}_${id}_raw.bw") into (bw4trans_ch, bw4seacr_ch, bw4macs2_ch)

        script:
        """
        samtools index ${sort_bam}
        bamCoverage -b ${sort_bam} -o ${group}_${id}_raw.bw
        """  
    }
}

/*
=======================================================================================
    Assess replicate reproducibility
=======================================================================================
*/

process assess_reproducibility {
    publishDir "${params.result}/alignment/bed", mode: "copy"

    input:
    tuple val(group), val(cond), val(id), path(frag_bed) from fragBed_assess_ch

    output:
    path("${group}_${id}_*.bed") into fragcount_bed_ch

    script:
    """
    awk ${params.args.awk_binLen} ${frag_bed} | ${params.args.assess_reproducibility} > ${group}_${id}_fragmentsCount.bin${params.args.binLen}.bed
    """
}

seqdepth_4_R_ch
    .map{group, cond, id, seqdepth -> [seqdepth]}
    .set{seqDepth_file_ch}

/*
=======================================================================================
    Alignment summary report
=======================================================================================
*/

process align_summary_R {
    publishDir "${params.result}/report_R", mode: "copy"

    input:
    path(hg38) from bowtie2_report_hg38_ch.collect()
    path(spikeIn) from bowtie2_report_spikeIn_ch.collect()
    path(picard) from rmdup_report_ch.collect()
    path(fragLen) from fragmentLen_ch.collect()
    path(bed) from fragcount_bed_ch.collect()
    path(seqDepth) from seqDepth_file_ch.collect()
    
    output:
    path("Alignment_summary.csv") into (alignsummary4Seacr_ch, alignsummary4Macs2_ch)
    path("Sequencing_Depth_Summary.png")
    path("Duplicate_Summary.png")
    path("fragment_size_summary.png")
    path("Replicate_Reproducibility.png")
    path("Sequencing_depth.png")
    
    script:
    """
    Rscript "$projectDir/alignment_report.r" $hg38 $spikeIn $picard $fragLen $bed $seqDepth
    """
}

/*
=======================================================================================
    Peak Calling
=======================================================================================
*/

if (params.use_ctrl) {
    if (params.workflow.seacr) {
        bedgraph_ch
            .branch {group, cond, id, path ->
                ctrl: cond == "ctrl"
                exp: cond == "exp"
            }
            .set {branch_bedgraph_ch}

        branch_bedgraph_ch.ctrl
            .combine(branch_bedgraph_ch.exp, by: 0)
            .set {seacr_ch}

       process seacr_with_ctrl {
            publishDir "${params.result}/peakcalling/seacr", mode: "copy"

            input:
            tuple val(group), val(ctrl), val(ctrl_id), path(ctrl_file), val(exp), val(exp_id), path(exp_file) from seacr_ch 

            output:
            tuple val(group), val(exp_id), path("${group}_${exp_id}_seacr_ctrl_peaks.*") into (seacr_ctrl4hp_ch, seacr_ctrl4report_ch)
            tuple val(group), val(exp_id), path("${group}_${exp_id}_seacr_top${params.args.seacr_top}_peaks.*") into seacr_top4report_ch

            script:
            """
            bash ${params.seacr_sh} ${exp_file} ${ctrl_file} ${params.args.seacr_norm_mode} ${group}_${exp_id}_seacr_ctrl_peaks
            bash ${params.seacr_sh} ${exp_file} ${params.args.seacr_top} ${params.args.seacr_norm_mode} ${group}_${exp_id}_seacr_top${params.args.seacr_top}_peaks
            """
       }
        
        seacr_ctrl4report_ch
            .map{group, id, file -> [file]}
            .set{seacr_ctrl4r_ch}
        seacr_top4report_ch
            .map{group, id, file -> [file]}
            .set{seacr_top4r_ch}
        bam4seacr_report_ch
            .map{group,cond, id, bam -> [bam]}
            .set{seacr_bam_ch}
        
        process seacr_with_ctrl_R {
            publishDir "${params.result}/report_R", mode: 'copy'

            input:
            path(peak_ctrl) from seacr_ctrl4r_ch.collect()
            path(peak_top) from seacr_top4r_ch.collect()
            path(bam) from seacr_bam_ch.collect()
            path(alignSummary) from alignsummary4Seacr_ch

            output:

            path("seacr_peak_summary.csv")
            path("seacr_frip_summary.csv")
            path("seacr_frip_summary.png")
        
            script:
            """
            Rscript "$projectDir/peak_seacr.r" $peak_ctrl $peak_top $bam $alignSummary
            """
        }
    }

    if (params.workflow.macs2) {
        bam4macs2_ch
            .branch {group, cond, id, path ->
                ctrl: cond == "ctrl"
                exp: cond == "exp"
            }
            .set {branch_macs2_ch}
        branch_macs2_ch.ctrl
            .combine(branch_macs2_ch.exp, by:0)
            .set {macs2_ch}

        process macs2 {
            publishDir "${params.result}/peakcalling/macs2", mode: "copy"

            input:
            tuple val(group), val(ctr), val(ctr_id), path(ctrl_bam), val(exp), val(exp_id), path(exp_bam) from macs2_ch

            output:
            tuple val(group), val(exp_id), path("${group}_${exp_id}.macs2_*_summary.txt")
            tuple val(group), val(exp_id), path("${group}_${exp_id}_macs2_*_peaks.xls")
            tuple val(group), val(exp_id), path("${group}_${exp_id}_macs2_topq*_summits.bed")

            tuple val(group), val(exp_id), path("${group}_${exp_id}_macs2_baseq*_summits.bed") into macs2_4hp_ch
            tuple val(group), val(exp_id), path("${group}_${exp_id}_macs2_*_peaks.narrowPeak") into macs2_4report_ch

            script:
            """
            macs2 callpeak -t ${exp_bam} -c ${ctrl_bam} ${params.args.macs2_args} -q ${params.args.macs2_base_qvalue} -n ${group}_${exp_id}_macs2_baseq${params.args.macs2_base_qvalue} 2>${group}_${exp_id}.macs2_base_summary.txt
            macs2 callpeak -t ${exp_bam} -c ${ctrl_bam} ${params.args.macs2_args} -q ${params.args.macs2_top_qvalue} -n ${group}_${exp_id}_macs2_topq${params.args.macs2_top_qvalue} 2>${group}_${exp_id}.macs2_top_summary.txt

            """
        }

        macs2_4report_ch
            .map{group, id, file -> file}
            .set{macs2_4r_ch}
        
        bam4macs2_report_ch
            .map{group,cond, id, bam -> bam}
            .set{macs2_bam_ch}
        
        process macs2_R {
            publishDir "${params.result}/report_R", mode: 'copy'

            input:
            path(peak) from macs2_4r_ch.collect()
            path(bam) from macs2_bam_ch.collect()
            path(alignSummary) from alignsummary4Macs2_ch

            output:

            path("macs2_peak_summary.csv")
            path("macs2_frip_summary.csv")
            path("macs2_frip_summary.png")
        
            script:
            """
            Rscript "$projectDir/peak_macs2.r" $peak $bam $alignSummary
            """
        }
    }
} else {
    process seacr_without_ctrl {
        publishDir "${params.result}/peakcalling/seacr", mode: "copy"

        input:
        tuple val(group), val(cond), val(id), path(bedgraph) from bedgraph_ch

        output:
        tuple val(group), val(id), path("${group}_${id}_seacr_noctrl_top*_peaks.*") into (seacr_noctrl4hp_ch, seacr_noctrl4report_ch)
        tuple val(group), val(id), path ("${group}_${id}_seacr_top${params.args.seacr_top}_peaks.*") into seacr_top4report_ch

        script:
        """
        bash ${params.seacr_sh} ${bedgraph} ${params.args.seacr_threshold_noconrl} ${params.args.seacr_norm_mode} ${group}_${id}_seacr_noctrl_top${params.args.seacr_threshold_noconrl}_peaks
        bash ${params.seacr_sh} ${bedgraph} ${params.args.seacr_top} ${params.args.seacr_norm_mode} ${group}_${id}_seacr_top${params.args.seacr_top}_peaks
        """
    }

    seacr_noctrl4report_ch
        .map{group, id, file -> [file]}
        .set{seacr_noctrl4r_ch}
    seacr_top4report_ch
        .map{group, id, file -> [file]}
        .set{seacr_top4r_ch}
    bam4seacr_report_ch
        .map{group,cond, id, bam -> [bam]}
        .set{seacr_bam_ch}

    process seacr_without_ctrl_report {
    publishDir "${params.result}/report_R", mode: 'copy'

    input:
    path(peak_noctrl) from seacr_noctrl4r_ch.collect()
    path(peak_top) from seacr_top4r_ch.collect()
    path(bam) from seacr_bam_ch.collect()
    path(alignSummary) from alignsummary4Seacr_ch

    output:
    path("seacr_peak_summary.csv")
    path("seacr_frip_summary.csv")
    path("seacr_frip_summary.png")

    script:
    """
    Rscript "$projectDir/peak_seacr.r" $peak_noctrl $peak_top $bam $alignSummary
    """
    } 
}

/*
=======================================================================================
    Visualization
=======================================================================================
*/

if (params.use_ctrl) {
    if (params.workflow.seacr) {
        seacr_4_summit_ch = seacr_ctrl4hp_ch
    } else { seacr_4_summit_ch = Channel.empty()}
} else {seacr_4_summit_ch = seacr_noctrl4hp_ch}


if (params.workflow.heatmap_on_trans) {
    Channel.fromPath(params.gene_for_heatmap).set{ref_heatmap_ch}
    
    bw4trans_ch
        .map {group, cond, id, path -> [path]}
        .collect()
        .set {bw4trans_collect_ch}

    process trans_heatmap {
        publishDir "${params.result}/heatmap/matrix", pattern:"*.mat.gz", mode:"copy"
        publishDir "${params.result}/heatmap", pattern:"*.png", mode:"copy"
        
        input:
        path(ref) from ref_heatmap_ch
        path(bw) from bw4trans_collect_ch
        
        output:
        path("matrix_gene.mat.gz")
        path("Over_Transcription_Units.png")
        
        script:
        """
        computeMatrix scale-regions -p ${params.threads} -S ${bw} -R ${ref} -o matrix_gene.mat.gz ${params.args.matrix_transcription}
        plotHeatmap -m matrix_gene.mat.gz -out Over_Transcription_Units.png ${params.args.plot_transcription}
        """
    }
}

if (params.workflow.heatmap_on_peak) {
    if (params.use_ctrl) {
        if (params.workflow.seacr) {
            
            process seacr_summitregion{
                publishDir "${params.result}/peakcalling/seacr", mode:"copy"

                input:
                tuple val(group), val(id), path(ctrl_peak) from seacr_ctrl4hp_ch

                output:
                tuple val(group), val(id), path("${group}_${id}_seacr_control_peaks.summitRegion.bed") into seacr_summit_ch

                script:
                """
                awk ${params.args.seacr_summit_awk} ${ctrl_peak} > ${group}_${id}_seacr_control_peaks.summitRegion.bed
                """
            }

            seacr_summit_ch
                .join(bw4seacr_ch, by:[0,1])
                .set {bw_seacr_ch}

            process seacr_withCtrl_heatmap {
                publishDir "${params.result}/heatmap/matrix", pattern:"*.mat.gz", mode:"copy"
                publishDir "${params.result}/heatmap", pattern:"*.png", mode:"copy"

                input:
                tuple val(group), val(id), path(sr_bed), path(bw) from bw_seacr_ch

                output:
                tuple val(group), val(id), path("${group}_${id}_seacr_withCtrl.mat.gz")
                tuple val(group), val(id), path("${group}_${id}_seacr_withCtrl.png")

                script:
                """
                computeMatrix reference-point -p ${params.threads} -S ${bw} -R ${sr_bed} -o ${group}_${id}_seacr_withCtrl.mat.gz ${params.args.matrix_peak}
                plotHeatmap -m ${group}_${id}_seacr_withCtrl.mat.gz -out ${group}_${id}_seacr_withCtrl.png ${params.args.plot_peak} --samplesLabel "${id} group: ${group}"
                """
            }
        }

        if (params.workflow.macs2) {

            macs2_4hp_ch
                .join(bw4macs2_ch, by:[0,1])
                .set {bw_macs2_ch}

            process macs2_heatmap {
                publishDir "${params.result}/heatmap/matrix", pattern:"*_macs2.mat.gz", mode:"copy"
                publishDir "${params.result}/heatmap", pattern:"*_macs2.png", mode:"copy"

                input:
                tuple val(group), val(id), path(macs_bed), path(bw) from bw_macs2_ch

                output:
                tuple val(group), val(id), path("${group}_${id}_macs2.mat.gz")
                tuple val(group), val(id), path("${group}_${id}_macs2.png")

                script:
                """
                computeMatrix reference-point -p ${params.threads} -S ${bw} -R ${macs_bed} -o ${group}_${id}_macs2.mat.gz ${params.args.matrix_peak}
                plotHeatmap -m ${group}_${id}_macs2.mat.gz -out ${group}_${id}_macs2.png ${params.args.plot_peak} --samplesLabel "${id} group: ${group}"
                """
            }
        }
    } else {
        process seacr_noCtrl_summitregion{
            publishDir "${params.result}/peakcalling/seacr", mode:"copy"

            input:
            tuple val(group), val(id), path(ctrl_peak) from seacr_noctrl4hp_ch

            output:
            tuple val(group), val(id), path("${group}_${id}_seacr_noctrl_peaks.summitRegion.bed") into seacr_summit_ch

            script:
            """
            awk ${params.args.seacr_summit_awk} ${ctrl_peak} > ${group}_${id}_seacr_nortcl_peaks.summitRegion.bed
            """
        }

        seacr_summit_ch
            .join(bw4seacr_ch, by:[0,1])
            .set {bw_seacr_ch}

        process seacr_noCtrl_heatmap {
            publishDir "${params.result}/heatmap/matrix", pattern:"*_seacr_noCtrl.mat.gz", mode:"copy"
            publishDir "${params.result}/heatmap", pattern:"*_seacr_noCtrl.png", mode:"copy"

            input:
            tuple val(group), val(id), path(sr_bed), path(bw) from bw_seacr_ch

            output:
            tuple val(group), val(id), path("${group}_${id}_seacr_noCtrl.mat.gz")
            tuple val(group), val(id), path("${group}_${id}_seacr_noCtrl.png")

            script:
            """
            computeMatrix reference-point -p ${params.threads} -S ${bw} -R ${sr_bed} -o ${group}_${id}_seacr_noCtrl.mat.gz ${params.args.matrix_peak}
            plotHeatmap -m ${group}_${id}_seacr_noCtrl.mat.gz -out ${group}_${id}_seacr_noCtrl.png ${params.args.plot_peak} --samplesLabel "${id} group: ${group}"
            """
        }
    }
}

