# Workflow for preprocessing BAM files for variant calling
from os.path import join

# TODO: Combine XX and XY to each rule. Create GATK within the environment.

configfile: "TCGA_LIHC_v2.config.json"

# Directories
SORTED_BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/TCGA_LIHC/WGS/GRCh38_sorted_BAM/" # path to directory for sorted BAM alignment files

# Samples
SAMPLES = config["TCGA_LIHC_WGS_without_rerun_samples"] #** Changed from TCGA_LIHC_WGS_samples
XX_SAMPLES = config["TCGA_LIHC_WGS_without_rerun_females"] #** Changed from, TCGA_LIHC_WGS_females
XY_SAMPLES = config["TCGA_LIHC_WGS_males"]

# Reference files
REF = config["XY_wholegenome_GRCh38_ref_with_viral_genomes_path"]
DBSNP = config["dbsnp_vcf"]
INTERVAL_LIST_XY = config["XY_interval_list"]
INTERVAL_LIST_XX = config["XX_interval_list"]


ruleorder: picard_add_readgroups > index_rg_bam > picard_mark_duplicates > index_picard_bam > gatk_baseq_recalibration > gatk_create_recalibrated_bam > index_processed_bam

rule all:
	input:
		# Defining the files that snakemake will attempt to produce as an output.
		# If there is no rule defined to produce the file, or if the file already
		# exists, snakemake will throw "Nothing to be done"
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XY_withoytYpar_HISAT2_aligned_sortedbycoord.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_recal_data.table", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_recal_data.table", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES)


rule picard_add_readgroups:
	input:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoytYpar_HISAT2_aligned_sortedbycoord.bam")
	output:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs.bam")
	params:
		SAMPLENAME = "{sample}",
	message: "Adding readgroups to {input.BAM_XX} and {input.BAM_XY}."
	run:
		if input.BAM_XX:
			shell("picard AddOrReplaceReadGroups I={input.BAM_XX} O={output.BAM_XX} RGID=0 RGLB=lib1 RGPL=Illumina RGPU=unit1 RGSM={params.SAMPLENAME}")
		if input.BAM_XY:
			shell("picard AddOrReplaceReadGroups I={input.BAM_XY} O={output.BAM_XY} RGID=0 RGLB=lib1 RGPL=Illumina RGPU=unit1 RGSM={params.SAMPLENAME}")


rule index_rg_bam:
	input:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs.bam")
	output:
		os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs.bam.bai"), #Missing from rule all
		os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs.bam.bai")
	message:
		"Indexing BAM file {input.BAM_XX} and {input.BAM_XY} with Samtools."
	params:
	run:
		for x in input:
			shell("samtools index {x}")


rule picard_mark_duplicates:
	input:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs.bam")
	output:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam"),
		METRICS = os.path.join(SORTED_BAM_AL_DIR,  "{sample}_marked_dup_metrics.txt") #Missing from rule all
	params:
		SAMPLENAME = "{sample}",
	message: "Removing duplicates from {input.BAM_XX} and {input.BAM_XY}."
	run:
		if input.BAM_XX:
			shell("picard MarkDuplicates I={input.BAM_XX} O={output.BAM_XX} M={output.METRICS} REMOVE_SEQUENCING_DUPLICATES=true")
		if input.BAM_XY:
			shell("picard MarkDuplicates I={input.BAM_XY} O={output.BAM_XY} M={output.METRICS} REMOVE_SEQUENCING_DUPLICATES=true")


rule index_picard_bam:
	input:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam")
	output:
		os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam.bai"), 
		os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam.bai")
	message:
		"Indexing BAM file {input.BAM_XX} and {input.BAM_XY} with Samtools."
	params:
	run:
		for x in input:
			shell("samtools index {x}")


rule gatk_baseq_recalibration:
	input:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam"),
		REF = REF,
		DBSNP = DBSNP,
		INTERVALS_XX = INTERVAL_LIST_XX,
		INTERVALS_XY = INTERVAL_LIST_XY,
		BAI_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam.bai"),
		BAI_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam.bai")
	output:
		DATA = os.path.join(SORTED_BAM_AL_DIR, "{sample}_recal_data.table")
	message: "Recalibrating base quality scores in {input.BAM_XX} and {input.BAM_XY} with GATK BaseRecalibrator"
	params:
		xmx = "16g"
	run:
		if input.BAM_XX:
			shell(
				"""
				gatk-launch BaseRecalibrator\
				--java-options "-Xmx{params.xmx}" \
				--reference {input.REF} \
				--input {input.BAM_XX} \
				--intervals {input.INTERVALS_XX} \
				--known-sites {input.DBSNP} \
				--output {output.DATA}
				""")
		if input.BAM_XY:
			shell(
				"""
				gatk-launch BaseRecalibrator\
				--java-options "-Xmx{params.xmx}" \
				--reference {input.REF} \
				--input {input.BAM_XY} \
				--intervals {input.INTERVALS_XY} \
				--known-sites {input.DBSNP} \
				--output {output.DATA}
				""")


rule gatk_create_recalibrated_bam:
	input:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups.bam"),
		BASEQ_TABLE = os.path.join(SORTED_BAM_AL_DIR, "{sample}_recal_data.table"),
		REF = REF,
		INTERVALS_XX = INTERVAL_LIST_XX,
		INTERVALS_XY = INTERVAL_LIST_XY
	output:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam")
	message: "Applying base quality scores in {input.BAM_XX} and {input.BAM_XY} with GATK ApplyBQSR"
	params:
		xmx = "16g"
	run:
		if input.BAM_XX:
			shell(
				"""
				gatk-launch ApplyBQSR\
				--java-options "-Xmx{params.xmx}" \
				--reference {input.REF} \
				--input {input.BAM_XX} \
				--bqsr-recal-file {input.BASEQ_TABLE} \
				--intervals {input.INTERVALS_XX} \
				--output {output.BAM_XX}
				""")
		if input.BAM_XY:
			shell(
				"""
				gatk-launch ApplyBQSR\
				--java-options "-Xmx{params.xmx}" \
				--reference {input.REF} \
				--input {input.BAM_XY} \
				--bqsr-recal-file {input.BASEQ_TABLE} \
				--intervals {input.INTERVALS_XY} \
				--output {output.BAM_XY}
				""")
	

rule index_processed_bam:
	input:
		BAM_XX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam"),
		BAM_XY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam")
	output:
		os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam.bai"),
		os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoutYpar_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam.bai")
	message: "Indexing BAM file {input.BAM_XX} and {input.BAM_XY} with Samtools."
	run:
		for x in input:
			shell("samtools index {x}")