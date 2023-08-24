#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process queryNCBI {
	tag "$run"
	label 'sratoolkit'

	input:
	val(run)
	val(results)
	val flag_1
	val flag_2

	output:
	tuple val(run), file("*.fastq.gz")

	when:
	((!mf.checkFile("$results/FASTQ/RAW", run, "q.gz") && \
	!mf.checkFile("$results/FASTQ/FILERED", run, "q.gz") && \
	!mf.checkFile("$results/BAM/", run, "bam") && \
	!mf.checkFile("$results/VCF/FILTERED", run, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", run, "vcf.gz")) && !params.help && !params.dry) || mf.checkFORCE('DOWNLOAD', params.FORCE)

	script:
	"""
	prefetch $run -N 1000

	fasterq-dump $run --skip-technical
	pigz *.fastq --best

	rm -rf $results/FASTQ/RAW/${run}*.fastq.gz
	mv *.fastq.gz $results/FASTQ/RAW
	ln -s $results/FASTQ/RAW/$run* .

	rm -rf $run
	"""
}

process fasta2NGS {
	tag "$sample"
	label 'samtools'

	input:
	tuple val(sample), val(compression), file(fasta)
	val(results)
	val flag_1
	val flag_2

	output:
	tuple val(sample), file("${sample}_1.fastq.gz"), file("${sample}_2.fastq.gz")

	when:
	((!mf.checkFile("$results/FASTQ/RAW", run, "q.gz") && \
	!mf.checkFile("$results/FASTQ/FILTERED", run, "q.gz") && \
	!mf.checkFile("$results/BAM/", run, "bam") && \
	!mf.checkFile("$results/VCF_FILTERED", run, "vcf.gz") && \
	!mf.checkFile("$results/VCF_RAW", run, "vcf.gz")) && !params.help && !params.dry) || mf.checkFORCE('DOWNLOAD', params.FORCE)

	script:
	if(compression == true)
		"""
		#generate 1m reads, 250 bp, do not add mutations, haploid
		zcat $fasta > unzipped.fasta
		wgsim unzipped.fasta ${sample}_1.fastq ${sample}_2.fastq -N 1000000 -1 250 -2 250 -h 1 -e 0 -r 0 -R 0 -X 0 -S 11031925
		bgzip ${sample}_1.fastq
		bgzip ${sample}_2.fastq

		rm -rf $results/FASTQ/RAW/${sample}*.fastq.gz
		mv ${sample}_1.fastq.gz $results/FASTQ/RAW ; mv ${sample}_2.fastq.gz $results/FASTQ/RAW
		ln -s $results/FASTQ/RAW/${sample}_1.fastq.gz . ; ln -s $results/FASTQ/RAW/${sample}_1.fastq.gz .
		"""
	else if(compression == false)
		"""
		#pretty much the same
		wgsim $fasta ${sample}_1.fastq ${sample}_2.fastq -N 1000000 -1 250 -2 250 -h 1 -e 0 -r 0 -R 0 -X 0 -S 11031925
		bgzip ${sample}_1.fastq
		bgzip ${sample}_2.fastq

		rm -rf $results/FASTQ/RAW/${sample}*.fastq.gz
		mv ${sample}_1.fastq.gz $results/FASTQ/RAW ; mv ${sample}_2.fastq.gz $results/FASTQ/RAW
		ln -s $results/FASTQ/RAW/${sample}_1.fastq.gz . ; ln -s $results/FASTQ/RAW/${sample}_1.fastq.gz .
		"""
	else
		error "You should not be there !"

}

mf = new myFunctions()

workflow download {

	results = file(params.results)
	file("$results/FASTQ/RAW").mkdirs()

	main:
	//import user fastq
	old_fastq = Channel.fromPath(params.input+"/*.{fq,fastq}.gz", followLinks: true)
		.map{it -> [it.simpleName.split("_1\$|_2\$|_R1\$|_R2\$")[0], it]}.groupTuple()
		.branch{
			paired: it[1].size() == 2
			single: it[1].size() == 1
		}
	old_fastq.paired.mix(old_fastq.single).map{it -> it[1]}.flatten().subscribe{ it -> mf.createSymLink(it.toString(), results.toString()+"/FASTQ/RAW") }

	old_single = old_fastq.single.map{it -> [ it[0], it[1][0] ] }
	old_paired = old_fastq.paired.map{it -> [ it[0], it[1][0], it[1][1] ] }

	//split the SRA list and channel it
	SRA = Channel.from(file(params.sra)).splitText().map{it -> it.trim()}.filter( it -> it !=~ /^#/ ) //when nohelp nodry

	//query NCBI and download the fastq
	queryNCBI(SRA, results, old_paired.collect().ifEmpty(true), old_single.collect().ifEmpty(true))

	queryNCBI.out.branch{
		paired: it[1].size() == 2
		single: it[1].size() == 0
	}.set{ new_fastq }

	new_single = new_fastq.single.map{ it -> [ it[0], it[1] ] }
	new_paired = new_fastq.paired.map{ it -> [ it[0], it[1][0], it[1][1] ] }

	//generate artificial fastq from fasta file
	fasta = Channel.fromPath(params.input+"/*.{fna,fa,fasta}.gz", followLinks: true).map{it -> [ it.simpleName, true, it]}
			.mix(Channel.fromPath(params.input+"/*.{fna,fa,fasta}", followLinks: true).map{it -> [ it.simpleName, false, it]})

	fasta2NGS(fasta, results, old_paired.collect().ifEmpty(true), old_single.collect().ifEmpty(true))
	fasta2NGS.out.map{it -> [ it[1], it[2] ] }.flatten().subscribe{ it -> mf.createSymLink(it.toString(), results.toString()+"/FASTQ/RAW") }

	//collision if downloaded is also a fasta but risk is minimal
	old_single.mix(old_single, old_paired, new_single, new_paired, fasta2NGS.out)
		.branch{
			paired: it.size() == 3
			single: it.size() == 2
		}.set{ fastq }

	emit:
	all_single_fastq = fastq.single
	all_paired_fastq = fastq.paired
}
