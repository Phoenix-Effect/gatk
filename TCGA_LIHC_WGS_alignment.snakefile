#! importing join
from os.path import join

# Workflow for aligning DNA sequence data to a sex specific reference with
# HISAT2, sorting and indexing BAM files with Samtools.
configfile: "TCGA_LIHC_v1.11.config.json"

# Tools
HISAT2 = "hisat2"
SAMTOOLS = "samtools"

# Reference genome files: XX with Y chromosome masked, XY with both Y-chromosomal PAR masked
# XX_HISAT2_INDEX_WITH_VIRAL_REF = config["XX_HISAT2_index_GRCh38_ref_with_viral_genomes"]
# XY_WITHOUTYPAR_HISAT2_INDEX_WITH_VIRAL_REF = config["XY_withoutYpar_HISAT2_index_GRCh38_ref_with_viral_genomes"]
XX_HISAT2_INDEX_WITH_VIRAL_REF = config["XX_GRCh38_ref_with_viral_genomes_HISAT2_index"]
XY_WITHOUTYPAR_HISAT2_INDEX_WITH_VIRAL_REF = config["XY_withoutYpar_GRCh38_ref_with_viral_genomes_HISAT2_index"]

# Directories
FQ_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/TCGA_LIHC/WGS/trimmed_fastqs_leading25_trailing25_winqual25//" # path to directory with trimmed FASTQ files
SAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/TCGA_LIHC/WGS/GRCh38_SAM/" # path to directory for SAM alignment files
BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/TCGA_LIHC/WGS/GRCh38_BAM/" # path to directory for BAM alignment files
SORTED_BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/TCGA_LIHC/WGS/GRCh38_sorted_BAM/" # path to directory for sorted BAM alignment files

# Samples
XX_SAMPLES = config["TCGA_LIHC_WGS_without_rerun_females"]
XY_SAMPLES= config["TCGA_LIHC_WGS_males"]
SAMPLES = config["TCGA_LIHC_WGS_without_rerun_samples"]

rule all:
	input:
		# Defining the files that snakemake will attempt to produce as an output.
		# If there is no rule defined to produce the file, or if the file already
		# exists, snakemake will throw "Nothing to be done"
		expand(SAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned.sam", SAM_AL_DIR=SAM_AL_DIR, sample=XX_SAMPLES),
		expand(BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned.bam", BAM_AL_DIR=BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SAM_AL_DIR + "{sample}_WGS_XY_withoytYpar_HISAT2_aligned.sam", SAM_AL_DIR=SAM_AL_DIR, sample=XY_SAMPLES),
		expand(BAM_AL_DIR + "{sample}_WGS_XY_withoytYpar_HISAT2_aligned.bam", BAM_AL_DIR=BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XY_withoytYpar_HISAT2_aligned_sortedbycoord.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_WGS_XY_withoytYpar_HISAT2_aligned_sortedbycoord.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES)


# TODO: use path.join
# Done
rule hisat2_xx_align_reads:
	input:
		R1 = os.path.join(FQ_DIR, "{sample}_trimmomatic_trimmed_paired_1.fastq"),
		R2 = os.path.join(FQ_DIR, "{sample}_trimmomatic_trimmed_paired_2.fastq")
	output:
		XX_SAM = os.path.join(SAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned.sam"),
		XY_SAM = os.path.join(SAM_AL_DIR, "{sample}_WGS_XY_withoytYpar_HISAT2_aligned.sam")
	params:
		XX_hisat2_index = XX_HISAT2_INDEX_WITH_VIRAL_REF,
		XY_hisat2_index = XY_WITHOUTYPAR_HISAT2_INDEX_WITH_VIRAL_REF,
		threads = 8
	message: "Mapping {wildcards.sample} reads to {params.XX_hisat2_index} and {params.XY_hisat2_index} with HISAT2."
	run:
		if input.R1:
			shell("hisat2 -q --phred33 -p {params.threads} -x {params.XX_hisat2_index} -s no -1 {input.R1} -2 {input.R2} -S {output.XX_SAM}")
		if input.R2:
			shell("hisat2 -q --phred33 -p {params.threads} -x {params.XY_hisat2_index} -s no -1 {input.R1} -2 {input.R2} -S {output.XY_SAM}")

rule xx_xy_sam_to_bam:
	input:
		SAMXX = os.path.join(SAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned.sam"),
		SAMXY = os.path.join(SAM_AL_DIR, "{sample}_WGS_XY_withoytYpar_HISAT2_aligned.sam")
	output:
		BAMXX = os.path.join(BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned.bam"),
		BAMXY = os.path.join(BAM_AL_DIR, "{sample}_WGS_XY_withoytYpar_HISAT2_aligned.bam")
	params:
	message: "Converting {input.SAMXX} and {input.SAMXY} to BAM, only outputting mapped reads."
	run:
		if input.SAMXX:
			shell("samtools view -b -F 4 {input.SAMXX} > {output.BAMXX}")
		if input.SAMXY:
			shell("samtools view -b -F 4 {input.SAMXY} > {output.BAMXY}")


rule xx_xy_sort_bam:
	input:
		BAMXX = os.path.join(BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned.bam"),
		BAMXY = os.path.join(BAM_AL_DIR, "{sample}_WGS_XY_withoytYpar_HISAT2_aligned.bam")
	output:
		SORTED_BAMXX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord.bam"),
		SORTED_BAMXY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoytYpar_HISAT2_aligned_sortedbycoord.bam")
	params:
	message: "Sorting BAM file {input.BAMXX} and {input.BAMXY}"
	run:
		if input.BAMXX:
			shell("samtools sort -O bam -o {output.SORTED_BAMXX} {input.BAMXX}")
		if input.BAMXY:
			shell("samtools sort -O bam -o {output.SORTED_BAMXY} {input.BAMXY}")

rule xx_xy_index_bam:
	input:
		BAMXX = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord.bam"),
		BAMXY = os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoytYpar_HISAT2_aligned_sortedbycoord.bam")
	output: 
		os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XX_HISAT2_aligned_sortedbycoord.bam.bai"),
		os.path.join(SORTED_BAM_AL_DIR, "{sample}_WGS_XY_withoytYpar_HISAT2_aligned_sortedbycoord.bam.bai")
	message: "Indexing BAM file {input.BAMXX} and {input.BAMXY} with Samtools."
	params:
	run:
		for x in input:
			shell("samtools index {x}")