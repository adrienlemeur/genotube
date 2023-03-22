#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process removeFastq {
	label 'script'

	input:
	tuple val(sample), file("${sample}.bam")
	val(results)

	when:
	!params.help && mf.checkFORCE('FASTQ', params.REMOVE)

	script:
	"""
	rm ${results}/FASTQ/RAW/${sample}* ${results}/FASTQ/TRIMMED/${sample}*
	"""
}

process removeBam {
	label 'script'

	input:
	tuple val(sample), file(vcf)
	val(results)

	when:
	!params.help && mf.checkFORCE('BAM', params.REMOVE)

	script:
	"""
	rm ${results}/BAM/RAW/${sample}*
	if test -f '$results/BAM/FILTERED/*.bam'; then
		rm ${results}/BAM/FILTERED/${sample}*
	fi
	"""
}

process removeVCF {
	label 'script'

	input:
	tuple val(sample), file(vcf), file("all_strain_info.txt")
	val(results)

	when:
	!params.help && mf.checkFORCE('VCF', params.REMOVE)

	script:
	"""
	rm ${results}/VCF/RAW/${sample}* ${results}/VCF/FILTERED/${sample}*
	"""
}

mf = new myFunctions()

workflow cleaner {
	take: all_mapping
	take: all_ann_vcf
	take: strain_info

	main:

	results = file(params.results)

	removeFastq(all_mapping, results)

	removeBam(all_ann_vcf, results)

	removeVCF(strain_info, results)

}