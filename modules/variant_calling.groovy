#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GATK {
	tag "$sample"
	label 'GATK'

	input:
	tuple val(sample), file(bam), file(bai)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}_GATK.vcf.gz")

	when:
	params.variant_caller == 'gatk' && !params.target_region && \
	(!mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") || mf.checkFORCE('CALL', params.FORCE) )

	script:
	"""
	gatk HaplotypeCaller \\
		--sample-ploidy 2 \\
		-R $fasta \\
		-I $bam -O ${sample}_GATK.vcf.gz \\
		--max-mnp-distance 1 \\
		--create-output-variant-index \\
		--native-pair-hmm-threads ${task.cpus}
	"""
}

process mpileup {
	tag "$sample"
	label 'bcftools'

	input:
	tuple val(sample), file(bam), file(bai)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}_samtools.vcf.gz")

	when:
	params.variant_caller == 'samtools' && !params.target_region && \
	(! mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") && \
	! mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") || mf.checkFORCE('CALL', params.FORCE) )

	script:
	"""
	bcftools mpileup -f $fasta $bam -Q 30 | \
		bcftools call --ploidy 2 -A -mv -Oz -o ${sample}_samtools.vcf.gz --threads ${task.cpus}
	"""
}

process freebayes {
	tag "$sample"
	label 'freebayes'

	input:
	tuple val(sample), file(bam), file(bai)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}.vcf.gz")

	when:
	params.variant_caller == 'freebayes' && !params.target_region && \
	(!mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('CALL', params.FORCE)

	script:
	"""
	freebayes -f $fasta $bam --vcf ${sample}_freebayes.vcf \
		--ploidy 2 -C 4 -m 30 \
		--report-all-haplotype-alleles \
		--use-reference-allele

	bgzip ${sample}_freebayes.vcf

	rm -rf $results/VCF/RAW/${sample}.vcf.gz
	mv ${sample}_freebayes.vcf.gz $results/VCF/RAW/${sample}.vcf.gz ; ln -s $results/VCF/RAW/${sample}.vcf.gz
	"""
}

process indexVCF {
	tag "$sample"
	label 'bcftools'

	input:
	tuple val(sample), file(vcf)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}.normed.vcf.gz"), file("${sample}.normed.vcf.gz.csi")

	when:
	!mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") || mf.checkFORCE('CALL', params.FORCE)

	script:
	"""
	bcftools index $vcf
	bcftools norm -f ${fasta} $vcf -Oz -o ${sample}.normed.vcf.gz

	bcftools index ${sample}.normed.vcf.gz
	"""
}

process fasta2Elfasta {
	label 'elprep'

	input:
	tuple file(fasta), file(fai), file(dict)	

	output:
	file("${fasta.baseName}.elfasta")

	when:
	params.target_region

	script:
	"""
	elprep fasta-to-elfasta $fasta ${fasta.simpleName}.elfasta
	"""
}

process onePassBamProcess {
	tag "$sample"
	label 'elprep'

	input:
	tuple val(sample), file(bam)
	file(elfasta)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}.vcf.gz")

	when:
	params.target_region && \
	(! mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") && \
	! mf.checkFile("$results/VCF/RAW", sample, "vcf.gz")) || mf.checkFORCE('CALL', params.FORCE)

	script:
	if(params.target_region){region = file(params.target_region)}
	"""
	elprep filter $bam ${sample}.dedup.bam \\
		--replace-read-group 'ID:$sample LB:LIB PL:Illumina PU:UNIT SM:$sample' \\
		--sorting-order coordinate \\
		--mark-duplicates --mark-optical-duplicates $results/ALL_REPORTS/BAM/DEDUP/${sample}_deduplication_metrics.txt \\
		--reference $elfasta \\
		--target-regions $region \\
		--reference-confidence NONE \\
		--haplotypecaller ${sample}.vcf.gz

	rm -rf $results/VCF/RAW/${sample}.vcf.gz
	mv ${sample}.vcf.gz $results/VCF/RAW/${sample}.vcf.gz
	ln -s $results/VCF/RAW/${sample}.vcf.gz ${sample}.vcf.gz

	rm -rf $results/BAM/${sample}.dedup.bam
	mv ${sample}.dedup.bam $results/BAM/${sample}.dedup.bam
	ln -s $results/BAM/${sample}.dedup.bam ${sample}.dedup.bam
	"""
}

process getDeletionRegion {
	label 'delly'

	input:
	tuple val(sample), file(bam), file("*.bai")
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}_RD.bcf")

	when:
	params.deletion_region && !mf.checkFile("$results/DELETION_REGION", sample, "vcf.gz")

	script:
	"""
	samtools index -b $bam
	delly call --svtype DEL -g $fasta $bam -q 20 -s 15 -o ${sample}_RD.bcf

	rm -rf $results/${sample}_RD.bcf

	mv ${sample}_RD.bcf $results/${sample}_RD.bcf
	ln -s $results/${sample}_RD.bcf ${sample}_RD.bcf
	"""
}


process bcfToVcf {
	label 'bcftools'

	input:
	tuple val(sample), file(bcf)
	val(results)

	output:
	tuple val(sample), file("${sample}_RD.vcf.gz")

	when:
	params.deletion_region && !mf.checkFile("$results/DELETION_REGION", sample, "vcf.gz")

	script:
	"""
	bcftools convert -Oz -o ${sample}_RD.vcf.gz $bcf
	rm -rf $results/${sample}_RD.vcf.gz
	ln ${sample}_RD.vcf.gz $results/${sample}_RD.vcf.gz
	rm $results/${sample}_RD.bcf
	"""
}

mf = new myFunctions()

workflow variant_calling {
	take: all_bam
	take: BAM_RAW
	take: index

	main:
		all_bam.view()

		results = file(params.results)
		file("$results/VCF/RAW").mkdirs()

		fasta2Elfasta(index)
		onePassBamProcess(BAM_RAW, fasta2Elfasta.out, index, results)

		GATK(all_bam, index, results)
		mpileup(all_bam, index, results)
		freebayes(all_bam, index, results)

		if( mf.checkFORCE('CALL', params.FORCE) ){
			old_vcf = Channel.empty()
		} else {
			old_vcf = Channel.fromPath(params.input+"/*.{vcf.gz,bcf}", followLinks: true)
				.map{it -> [ it.simpleName, it]}
			old_vcf.subscribe{ it -> mf.createSymLink(it[1].toString(), results.toString()+"/BAM") }

			old_vcf = Channel.fromPath(results+"/VCF/RAW/*.vcf.gz", followLinks: true)
				.map{it -> [ it.simpleName, it ] }
		}

		indexVCF(old_vcf.mix(freebayes.out, GATK.out, mpileup.out), index, results)

		all_vcf = indexVCF.out.mix(onePassBamProcess.out)

		//sort deletion regions in vcf files
		getDeletionRegion(all_bam, index, results+"/DELETION_REGION")
		bcfToVcf(getDeletionRegion.out, results+"/DELETION_REGION")

	emit:
		all_vcf
		garbage = indexVCF.out
}
