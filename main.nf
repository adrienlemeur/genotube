#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
------------------------------------------------------------------------------------------------------------
			 adrienlemeur/genotube-td
------------------------------------------------------------------------------------------------------------

	Genotube Analysis Pipeline. Started on 12-11-2021

	#### Homepage / Documentation
	https://github.com/adrienlemeur/genotube-td

	#### Authors
	Adrien Le Meur
	Fiona Hak
	with the help of Guislaine Refrégier & Rimma Zn
	#### Version : 4.0
	#### Name : Washed Ashore

------------------------------------------------------------------------------------------------------------
*/

// sub-workflow import
include { initialisation }	from		'./modules/initialisation.groovy'
include { download }		from		'./modules/download.groovy'
include { cleaner }		from		'./modules/cleaner.groovy'
include { process_fastq }	from		'./modules/process_fastq.groovy'
include { align }		from		'./modules/align.groovy'
include { index }		from		'./modules/indexing.groovy'
include { process_bam }	from			'./modules/process_bam.groovy'
include { variant_calling }	from		'./modules/variant_calling.groovy'
include { process_vcf }	from			'./modules/process_vcf.groovy'
include { profiling }		from		'./modules/profiling.groovy'
include { build_tree } 	from			'./modules/treebuild.groovy'
include { multiqc_report }	from		'./modules/multiqc_report.groovy'

mf = new myFunctions()

workflow {
	main:

		mf.checkParameters(params.variant_caller, params.contam_check, params.species, params.taxonomy, params.tree_build, params.FORCE, params.SKIP, params.REMOVE, params.download_K2_DB)

		initialisation()
		index()
		download()

		process_fastq(download.out, index.out.samtools_picard_index)

		align(process_fastq.out.all_single_trimmed, process_fastq.out.all_paired_trimmed, index.out.bwa_index, index.out.samtools_picard_index)
		process_bam(align.out.all_mapping, index.out.samtools_picard_index)

		variant_calling(process_bam.out.all_processed_bam, align.out.all_mapping, index.out.samtools_picard_index)

		process_vcf(variant_calling.out.all_vcf, process_bam.out.coverage_info, index.out.samtools_picard_index, index.out.snpeff_emit_signal)

		profiling(process_vcf.out.all_ann_vcf, process_vcf.out.all_raw_vcf, index.out.binExec_emit_signal, index.out.samtools_picard_index)

        	cleaner(align.out.all_mapping, process_vcf.out.all_ann_vcf, profiling.out.strain_info)

		build_tree(process_vcf.out.all_ann_vcf, profiling.out.strain_info, index.out.samtools_picard_index)
		multiqc_report(process_vcf.out.end_signal, process_vcf.out.all_ann_vcf)
}

