# Workflow for FastQC, MultiQC, and adapter trimming using Trimmomatic.
from os.path import join

configfile: "TCGA_LIHC_v1.11.config.json"

# Directory variables
fastq_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/TCGA_LIHC/WGS/stripped_fastqs/"
fastqc_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/TCGA_LIHC/WGS/fastqc/"
trimmed_fastqs = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/TCGA_LIHC/WGS/trimmed_fastqs_leading25_trailing25_winqual25/"
trimmed_fastqc_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/TCGA_LIHC/WGS/fastqc/trimmed_new/"
benchmark_dir = "/home/saghafoo/eqtl/bench/"

# Tools
fastqc_path = "fastqc"
multiqc_path = "multiqc"
trimmomatic_path = "trimmomatic"

SAMPLES = config["TCGA_LIHC_WGS_without_rerun_samples"]

ruleorder: fastqc_analysis > multiqc > trimmomatic > fastqc_analysis_trimmomatic_trimmed_paired > fastqc_analysis_trimmomatic_trimmed_paired > multiqc_trimmed_paired

rule all:
	input:
		# Defining the files that snakemake will attempt to produce as an output.
		# If there is no rule defined to produce the file, or if the file already
		# exists, snakemake will throw "Nothing to be done"
		expand(fastqc_directory + "{sample}_fq1_fastqc.html", fastqc_directory = fastqc_directory, sample=SAMPLES),
		expand(fastqc_directory + "{sample}_fq2_fastqc.html", fastqc_directory = fastqc_directory, sample=SAMPLES),
		# (fastqc_directory + "multiqc_report.html"),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq", trimmed_fastqs=trimmed_fastqs, sample=SAMPLES),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_1.fastq", trimmed_fastqs=trimmed_fastqs, sample=SAMPLES),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq", trimmed_fastqs=trimmed_fastqs, sample=SAMPLES),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_2.fastq", trimmed_fastqs=trimmed_fastqs, sample=SAMPLES),
		expand(trimmed_fastqs + "{sample}_trimmomatic.log", trimmed_fastqs=trimmed_fastqs, sample=SAMPLES),
		expand(trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_1_fastqc.html", sample=SAMPLES),
		expand(trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_2_fastqc.html", sample=SAMPLES),
		(trimmed_fastqc_directory + "multiqc_report.html")

# TODO: remove lambda, use path.join
# Done

rule fastqc_analysis:
	input:
		fq1 = os.path.join(fastq_directory, "{sample}_WGS_stripped_fq_1.fastq.gz"),
		fq2 = os.path.join(fastq_directory, "{sample}_WGS_stripped_fq_2.fastq.gz")
	output:
		fq1_zip =  os.path.join(fastqc_directory, "{sample}_fq1_fastqc.zip"),
		fq1_html = os.path.join(fastqc_directory, "{sample}_fq1_fastqc.html"),
		fq2_zip =  os.path.join(fastqc_directory, "{sample}_fq2_fastqc.zip"),
		fq2_html = os.path.join(fastqc_directory, "{sample}_fq2_fastqc.html")
	params:
		fastqc = fastqc_path,
		fastqc_dir = fastqc_directory,
		fq1_prefix = "{sample}_WGS_stripped_fq_1",
		fq2_prefix = "{sample}_WGS_stripped_fq_2"
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1};
		{params.fastqc} -o {params.fastqc_dir} {input.fq2};
		mv {params.fastqc_dir}{params.fq1_prefix}_fastqc.html {output.fq1_html};
		mv {params.fastqc_dir}{params.fq1_prefix}_fastqc.zip {output.fq1_zip};
		mv {params.fastqc_dir}{params.fq2_prefix}_fastqc.html {output.fq2_html};
		mv {params.fastqc_dir}{params.fq2_prefix}_fastqc.zip {output.fq2_zip}
		"""

rule multiqc:
	input:
	output:
		os.path.join(fastqc_directory, "multiqc_report.html")
	message: "Running MultiQC for FastQC reports located in {params.fastqc_dir}"
	params:
		fastqc_dir = fastqc_directory,
		output_dir = fastqc_directory
	shell:
		"""
		multiqc {params.fastqc_dir}*_fastqc.zip --outdir {params.output_dir} --interactive --verbose
		"""

# TODO: Change leading value
rule trimmomatic:
	input:
		fq1 = os.path.join(fastq_directory, "{sample}_WGS_stripped_fq_1.fastq.gz"),
		fq2 = os.path.join(fastq_directory, "{sample}_WGS_stripped_fq_2.fastq.gz"),
		ADAPTER_FASTA = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/adapter_sequences.fa"
	output:
		paired_1 =   os.path.join(trimmed_fastqs, "{sample}_trimmomatic_trimmed_paired_1.fastq"),
		unpaired_1 = os.path.join(trimmed_fastqs, "{sample}_trimmomatic_trimmed_unpaired_1.fastq"),
		paired_2 =   os.path.join(trimmed_fastqs, "{sample}_trimmomatic_trimmed_paired_2.fastq"),
		unpaired_2 = os.path.join(trimmed_fastqs, "{sample}_trimmomatic_trimmed_unpaired_2.fastq"),
		logfile =    os.path.join(trimmed_fastqs, "{sample}_trimmomatic.log")
	benchmark:
		"bench/{sample}.tsv"
	params:
		trimmomatic = trimmomatic_path,
		threads = 4,
		seed_mismatches = 2,
		palindrome_clip_threshold = 30,
		simple_clip_threshold = 10,
		leading = 25,
		trailing = 25,
		winsize = 4,
		winqual = 25,
		minlen = 80
	shell:
		"""
		{params.trimmomatic} PE -threads {params.threads} -phred33 -trimlog {output.logfile} \
		{input.fq1} {input.fq2} {output.paired_1} {output.unpaired_1} \
		{output.paired_2} {output.unpaired_2} \
		ILLUMINACLIP:{input.ADAPTER_FASTA}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} \
		LEADING:{params.leading} TRAILING:{params.trailing} \
		SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}
		"""


rule fastqc_analysis_trimmomatic_trimmed_paired:
	input:
		fq1 = os.path.join(trimmed_fastqs , "{sample}_trimmomatic_trimmed_paired_1.fastq"),
		fq2 = os.path.join(trimmed_fastqs , "{sample}_trimmomatic_trimmed_paired_2.fastq")
	output:
		html1 = os.path.join(trimmed_fastqc_directory, "{sample}_trimmomatic_trimmed_paired_1_fastqc.html"),
		html2 = os.path.join(trimmed_fastqc_directory, "{sample}_trimmomatic_trimmed_paired_2_fastqc.html")
	params:
		fastqc = fastqc_path,
		fastqc_dir = trimmed_fastqc_directory
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1} {input.fq2}
		"""


rule multiqc_trimmed_paired:
	input:
	output:
		os.path.join(trimmed_fastqc_directory, "multiqc_report.html"),
		os.path.join(trimmed_fastqc_directory, "multiqc_data")
	message: "Running MultiQC for post-trimming FastQC reports located in {params.fastqc_dir}"
	params:
		fastqc_dir = trimmed_fastqc_directory,
		output_dir = trimmed_fastqc_directory
	shell:
		"""
		multiqc {params.fastqc_dir}*trimmomatic*_fastqc.zip --outdir {params.output_dir} --interactive --verbose
		"""
