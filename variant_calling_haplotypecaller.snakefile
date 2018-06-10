# Workflow for calling germline SNP and indel variants using GATK HaplotypeCaller.
from os.path import join

configfile: "TCGA_LIHC_v2.config-singlefile.json"

# Reference files: reference genome sequences in FASTA format, index .fa.fai
# files created with Samtools faidx, sequence dictionary .dict files created
# with Picard CreateSequenceDictionary
REF = config["GRCh38_wholegenome_ref_path"]

# Database and interval files. Chromosome notation and names must match those
# in the input and reference files. As GATK does not currently support VCF
# version 4.2, VCF files must be formatted according to version 4.1.
DBSNP = config["dbsnp_vcf"] # Tabix indexed, bgzipped dbSNP VCF
XX_INTERVAL_LIST = config["chrXX_interval_list"]
XY_INTERVAL_LIST = config["chrXY_interval_list"]

# Directories
BAM_DIR = config["BAM_DIR"]
VCF_DIR = config["VCF_DIR"]
BENCHMARK_DIR = "/home/saghafoo/varcall/bench/"

# Samples
TCGA_LIHC_WES_normal_female_samples = config["TCGA_LIHC_WES_normal_female_samples"]
TCGA_LIHC_WES_normal_male_samples = config["TCGA_LIHC_WES_normal_male_samples"]

#ruleorder: genotype_maleXY_gvcfs_chr > genotype_male_gvcfs_chr

def list_to_variant_str(list):
	finalString = ''
	for i in list:
		finalString += '--variant ' + VCF_DIR + i + '_rawLikehoods.g.vcf '
	return finalString

def get_genotype_gvcfs_chr(type='all'):
	if type == 'female':
		return expand(VCF_DIR + "{sample}_rawLikehoods.g.vcf", sample=TCGA_LIHC_WES_normal_female_samples)
	if type == 'male':
		return expand(VCF_DIR + "{sample}_rawLikehoods.g.vcf", sample=TCGA_LIHC_WES_normal_male_samples)
	return expand(VCF_DIR + "{sample}_rawLikehoods.g.vcf", sample=TCGA_LIHC_WES_normal_male_samples + TCGA_LIHC_WES_normal_female_samples)


#	expand(VCF_DIR + "TCGA_LIHC_adjacent_male_genotype_gvcfs_chr{chr}.vcf", chr=list(map(str, (list(range(1,23)) + ['XandY'])))),
#	expand(VCF_DIR + "TCGA_LIHC_adjacent_female_genotype_gvcfs_chr{chr}.vcf", chr=list(map(str, (list(range(1,23)) + ['X']))))

rule all:
	input:
		get_genotype_gvcfs_chr('male')

rule gatk_haplotypecaller_gvcf:
	input:
		BAM = os.path.join(BAM_DIR, "{sample}_HISAT2_aligned_sortedbycoord_wRGs_rmDups_BQcalib.bam"),
		REF = REF,
		INTERVALSXY = XY_INTERVAL_LIST,
		INTERVALSXX = XX_INTERVAL_LIST,
		DBSNP = DBSNP
	output:
		VCF = os.path.join(VCF_DIR, "{sample}_rawLikehoods.g.vcf")
	log:
		"logs/gatk_haplotypecaller_gvcf_{sample}.log"
	params:
		xmx = "48g",
		xms = "48g",
		threads = "16",
		min_pruning = 2,
		heterozygosity=0.001, # Heterozygosity value used to compute prior likelihoods for any locus
		indel_heterozygosity=1.25E-4, # Heterozygosity for indel calling
		maxReadsInRegionPerSample=10000, # Maximum reads in an active region
		min_base_quality_score=10, # Minimum base quality required to consider a base for calling
		minReadsPerAlignmentStart=10, # Minimum number of reads sharing the same alignment start for each genomic location in an active region
		sample_ploidy=2, # Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)
		standard_min_confidence_threshold_for_calling=10.0, # The minimum phred-scaled confidence threshold at which variants should be called
		max_alternate_alleles=6, # Maximum number of alternate alleles to genotype
		output_mode="EMIT_ALL_CONFIDENT_SITES" # Which type of calls we should output: EMIT_VARIANTS_ONLY, EMIT_ALL_SITES
	message: "Calculating genotype likelihoods from {input.BAM} using GATK HaplotypeCaller and producing {output.VCF}"
	benchmark:
		BENCHMARK_DIR + "{sample}_bench.tsv"
	run:
		if "_XX_" in input.BAM:
			shell("""
				gatk-launch --java-options "-Xmx{params.xmx} -Xms{params.xms}" \
				HaplotypeCaller \
				--native-pair-hmm-threads {params.threads} \
				--emit-ref-confidence GVCF \
				--reference {input.REF} \
				--input {input.BAM} \
				--min-pruning {params.min_pruning} \
				--heterozygosity {params.heterozygosity} \
				--indel-heterozygosity {params.indel_heterozygosity} \
				--max-reads-per-alignment-start {params.maxReadsInRegionPerSample} \
				--min-base-quality-score {params.min_base_quality_score} \
				--sample-ploidy {params.sample_ploidy} \
				--standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold_for_calling} \
				--max-alternate-alleles {params.max_alternate_alleles} \
				--output-mode {params.output_mode} \
				--dbsnp {input.DBSNP} \
				--intervals {input.INTERVALSXX} \
				--output {output.VCF}
				""")
		if "_XY_" in input.BAM:
			shell("""
				gatk-launch --java-options "-Xmx{params.xmx} -Xms{params.xms}" \
				HaplotypeCaller \
				--native-pair-hmm-threads {params.threads} \
				--emit-ref-confidence GVCF \
				--reference {input.REF} \
				--input {input.BAM} \
				--min-pruning {params.min_pruning} \
				--heterozygosity {params.heterozygosity} \
				--indel-heterozygosity {params.indel_heterozygosity} \
				--max-reads-per-alignment-start {params.maxReadsInRegionPerSample} \
				--min-base-quality-score {params.min_base_quality_score} \
				--sample-ploidy {params.sample_ploidy} \
				--standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold_for_calling} \
				--max-alternate-alleles {params.max_alternate_alleles} \
				--output-mode {params.output_mode} \
				--dbsnp {input.DBSNP} \
				--intervals {input.INTERVALSXY} \
				--output {output.VCF}
				""")

# rule genotype_male_gvcfs_chr:
# 	input:
# 		REF = REF,
# 		INTERVALS = lambda wildcards: config['chr' + str(wildcards.chr) + '_interval_list'],
# 		CHRLIST = get_genotype_gvcfs_chr('male')
# 	output:
# 		VCF = os.path.join(VCF_DIR, "TCGA_LIHC_adjacent_male_genotype_gvcfs_chr{chr}.vcf") 
# 	params:
# 		GVCF_LIST_MALE = list_to_variant_str(TCGA_LIHC_WES_normal_male_samples), 
# 		xmx = "98g",
# 		xms = "98g",
# 		threads = "24",
# 		heterozygosity=0.001, # Heterozygosity value used to compute prior likelihoods for any locus
# 		indel_heterozygosity=1.25E-4, # Heterozygosity for indel calling
# 		sample_ploidy=2, # Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)
# 		standard_min_confidence_threshold_for_calling=10.0, # The minimum phred-scaled confidence threshold at which variants should be called
# 		max_alternate_alleles=6, # Maximum number of alternate alleles to genotype
# 	log:
# 		"logs/genotype_gvcfs_chr_{output.VCF}.log"
# 	message: "Merging gVCF records for {output.VCF}."
# 	run:
# 		shell(
# 			"""
# 			gatk-launch --java-options "-Xmx{params.xmx} -Xms{params.xms}" \
# 			GenotypeGVCFs \
# 			--reference {input.REF} \
# 			{params.GVCF_LIST_MALE} \
# 			--heterozygosity {params.heterozygosity} \
# 			--indel-heterozygosity {params.indel_heterozygosity} \
# 			--sample-ploidy {params.sample_ploidy} \
# 			--standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold_for_calling} \
# 			--max-alternate-alleles {params.max_alternate_alleles} \
# 			--intervals {input.INTERVALS} \
# 			--output {output.VCF}
# 			""")

# rule genotype_maleXY_gvcfs_chr:
# 	input:
# 		REF = REF,
# 		INTERVALSX = config['chrX_interval_list'],
# 		INTERVALSY = config['chrY_interval_list'],
# 		CHRLIST = get_genotype_gvcfs_chr('male')
# 	output:
# 		VCF = os.path.join(VCF_DIR, "TCGA_LIHC_adjacent_male_genotype_gvcfs_chrXandY.vcf") 
# 	params:
# 		GVCF_LIST_MALE = list_to_variant_str(TCGA_LIHC_WES_normal_male_samples), 
# 		xmx = "98g",
# 		xms = "98g",
# 		threads = "24",
# 		heterozygosity=0.001, # Heterozygosity value used to compute prior likelihoods for any locus
# 		indel_heterozygosity=1.25E-4, # Heterozygosity for indel calling
# 		sample_ploidy=2, # Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)
# 		standard_min_confidence_threshold_for_calling=10.0, # The minimum phred-scaled confidence threshold at which variants should be called
# 		max_alternate_alleles=6, # Maximum number of alternate alleles to genotype
# 	log:
# 		"logs/genotype_gvcfs_chr_{output.VCF}.log"
# 	message: "Merging gVCF records for {output.VCF}."
# 	run:
# 		shell(
# 			"""
# 			gatk-launch --java-options "-Xmx{params.xmx} -Xms{params.xms}" \
# 			GenotypeGVCFs \
# 			--reference {input.REF} \
# 			{params.GVCF_LIST_MALE} \
# 			--heterozygosity {params.heterozygosity} \
# 			--indel-heterozygosity {params.indel_heterozygosity} \
# 			--sample-ploidy {params.sample_ploidy} \
# 			--standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold_for_calling} \
# 			--max-alternate-alleles {params.max_alternate_alleles} \
# 			--intervals {input.INTERVALSX} \
# 			--intervals {input.INTERVALSY} \
# 			--output {output.VCF}
# 			""")


# rule genotype_female_gvcfs_chr:
# 	input:
# 		REF = REF,
# 		INTERVALS = lambda wildcards: config['chr' + str(wildcards.chr) + '_interval_list'],
# 		CHRLIST = get_genotype_gvcfs_chr('female')
# 	output:
# 		VCF = os.path.join(VCF_DIR, "TCGA_LIHC_adjacent_female_genotype_gvcfs_chr{chr}.vcf")
# 	params:
# 		GVCF_LIST_FEMALE = list_to_variant_str(TCGA_LIHC_WES_normal_female_samples),
# 		xmx = "98g",
# 		xms = "98g",
# 		threads = "24",
# 		heterozygosity=0.001, # Heterozygosity value used to compute prior likelihoods for any locus
# 		indel_heterozygosity=1.25E-4, # Heterozygosity for indel calling
# 		sample_ploidy=2, # Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)
# 		standard_min_confidence_threshold_for_calling=10.0, # The minimum phred-scaled confidence threshold at which variants should be called
# 		max_alternate_alleles=6, # Maximum number of alternate alleles to genotype
# 	log:
# 		"logs/genotype_gvcfs_chr_{output.VCF}.log"
# 	message: "Merging gVCF records for {output.VCF}."
# 	run:
# 		shell(
# 			"""
# 			gatk-launch --java-options "-Xmx{params.xmx} -Xms{params.xms}" \
# 			GenotypeGVCFs \
# 			--reference {input.REF} \
# 			{params.GVCF_LIST_FEMALE} \
# 			--heterozygosity {params.heterozygosity} \
# 			--indel-heterozygosity {params.indel_heterozygosity} \
# 			--sample-ploidy {params.sample_ploidy} \
# 			--standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold_for_calling} \
# 			--max-alternate-alleles {params.max_alternate_alleles} \
# 			--intervals {input.INTERVALS} \
# 			--output {output.VCF}
# 			""")