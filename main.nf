#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
------------------------------------------------------------------------------------------------------------
			 adrienlemeur/genotube-td
------------------------------------------------------------------------------------------------------------

	Genotube Analysis Pipeline. Started on 12-11-2021
	Last Update : 12-06-2023

	#### Homepage / Documentation
	https://github.com/adrienlemeur/genotube

	#### Authors
	Adrien Le Meur
	Fiona Hak
	with the help of Guislaine Refrégier & Rimma Zn
	#### Version : 6.0
	#### Name : Moses

------------------------------------------------------------------------------------------------------------
*/

// sub-workflow import
include { initialisation }	from		'./modules/initialisation.groovy'
include { download }		from		'./modules/download.groovy'
include { cleaner }		from		'./modules/cleaner.groovy'
include { process_fastq }	from		'./modules/process_fastq.groovy'
include { align }		from		'./modules/align.groovy'
include { index }		from		'./modules/indexing.groovy'
include { process_bam }		from		'./modules/process_bam.groovy'
include { variant_calling }	from		'./modules/variant_calling.groovy'
include { process_vcf }		from		'./modules/process_vcf.groovy'
include { profiling }		from		'./modules/profiling.groovy'
include { build_tree } 		from		'./modules/treebuild.groovy'
include { multiqc_report }	from		'./modules/multiqc_report.groovy'

mf = new myFunctions()

workflow {
	main:

		mf.checkParameters(params.variant_caller, params.contam_check, params.species, params.taxonomy, params.tree_build, params.FORCE, params.SKIP, params.REMOVE, params.download_K2_DB)

		initialisation()
		index()

		if(params.help || params.dry){
			exit(0)
		} else {

			//to do :
			//multiQC is still launched before the end of analysis
			//VCF and FILTERED VCF to do
			//VCF channel true error
			//TB-detective outputs only one file, switch back to several files or debug it
			//VCF can be duplicated (does not know if it still exists)
			//does not input bam check skip is not set

			download()
			process_fastq(download.out, index.out.samtools_picard)

			align(process_fastq.out.all_single_trimmed, process_fastq.out.all_paired_trimmed, index.out.bwa, index.out.samtools_picard)
			process_bam(align.out.all_mapping, index.out.samtools_picard)

			variant_calling(process_bam.out.all_processed_bam, align.out.all_mapping, index.out.samtools_picard)
			process_vcf(variant_calling.out.all_vcf, index.out.samtools_picard, index.out.snpeff_emit_signal)

			profiling(process_vcf.out.all_ann_vcf, index.out.binExec_emit_signal, index.out.samtools_picard)

			cleaner(align.out.all_mapping, process_vcf.out.all_ann_vcf, profiling.out.strain_info)

			build_tree(process_vcf.out.all_ann_vcf, profiling.out.strain_info, index.out.samtools_picard)
			multiqc_report(process_vcf.out.end_signal, process_vcf.out.all_ann_vcf)
		}
}

