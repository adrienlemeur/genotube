#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process spliceReference {
	//Create a fasta file for read filtering with bbduk
	//Note : the larger the target regions, the more memory consuming bbudk and elprep will be
	label 'bedtools'

	input:
	tuple file(fasta), file(fai), file(dict)

	output:
	tuple file("${fasta.simpleName}_spliced.fasta"), val(region.simpleName)

	when:
	params.target_region

	script:
	reference = fasta.simpleName
	if(params.target_region){region = file(params.target_region)}
	"""
	#get genome file for bedtools slop
	cut -f1,2 $fai > ${reference}.genome

	#add 200 nuc before and after every site for read filtering
	bedtools slop -i $region -b 300 -g ${reference}.genome > sloppy_bed.bed

	#copy reference header for new reference
	head -n1 $fasta > ${reference}_spliced.fasta

	#get a nucleotide sequence with only target regions +200/+200 nuc
	bedtools getfasta -fi $fasta -bed sloppy_bed.bed | grep -v ">" | tr -d '\n' >> ${reference}_spliced.fasta
	"""
}

process single_readCleaning {
	//bbduk kmer based filter
	tag "$sample"
	label 'bbtools'

	input:
	tuple val(sample), file(f1)
	tuple file(referenceSequence), val(region_name)
	val(results)

	output:
	tuple val("${sample}_${region_name}"), file("${sample}_${region_name}.filtered.fastq.gz")

	when:
	params.target_region && !params.dry && !params.help && (!mf.checkFile("$results/FASTQ/TRIMMED", sample, "q.gz") && \
	!mf.checkFile("$results/BAM/RAW", sample, "bam") && \
	!mf.checkFile("$results/BAM/FILTERED", sample, "bam") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('TRIM', params.FORCE)

	script:
	"""
	#kmer based filter
	bbduk.sh in=$f1 outm=f1.fastq.gz ref=$referenceSequence k=31 stats=${sample}_bbstats.txt -Xmx5g

	#filtering maybe unpair reads
	#repair faulty files
	repair.sh in=f1.fastq.gz out=${sample}_${region_name}.filtered.fastq.gz
	rm f1.fastq.gz

	#moves stats to the report folder for aggregation with multiqc
	mv ${sample}_bbstats.txt "$results/ALL_REPORTS/FASTQ/KMERFILTER/"
	"""
}

process paired_readCleaning {
	//bbduk kmer based filter
	tag "$sample"
	label 'bbtools'

	input:
	tuple val(sample), file(f1), file(f2)
	tuple file(referenceSequence), val(region_name)
	val(results)

	output:
	tuple val("${sample}_${region_name}"), file("${sample}_${region_name}_1.filtered.fastq.gz"), file("${sample}_${region_name}_2.filtered.fastq.gz")

	when:
	params.target_region && !params.dry && !params.help && (!mf.checkFile("$results/FASTQ/TRIMMED", sample, "q.gz") && \
	!mf.checkFile("$results/BAM/RAW", sample, "bam") && \
	!mf.checkFile("$results/BAM/FILTERED", sample, "bam") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('TRIM', params.FORCE)

	script:
	"""
	#kmer based filter
	bbduk.sh in=$f1 in2=$f2 outm=f1.fastq.gz outm2=f2.fastq.gz ref=$referenceSequence k=31 stats=${sample}_bbstats.txt -Xmx5g

	#filtering may unpair reads
	#repair faulty files
	repair.sh in=f1.fastq.gz in2=f2.fastq.gz out=${sample}_${region_name}_1.filtered.fastq.gz out2=${sample}_${region_name}_2.filtered.fastq.gz
	rm f1.fastq.gz f2.fastq.gz

	#moves stats to the report folder for aggregation with multiqc
	mv ${sample}_bbstats.txt "$results/ALL_REPORTS/FASTQ/KMERFILTER/"
	"""
}

process single_trimming {
	//Remove adaptaters, merge reads that can be merged
	tag "$sample"
	label 'fastp'

	input:
	tuple val(sample), file(f1)
	val(results)

	output:
	tuple val(sample), file("${sample}.trimmed.fastq.gz")

	when:
	!params.dry && !params.help && (!mf.checkFile("$results/FASTQ/TRIMMED", sample, "q.gz") && \
	!mf.checkFile("$results/BAM/RAW", sample, "bam") && \
	!mf.checkFile("$results/BAM/FILTERED", sample, "bam") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('TRIM', params.FORCE)

	script:
	"""
	fastp -i $f1 -o ${sample}.trimmed.fastq.gz \\
			-q 25 -l 20 -e 20 --dont_eval_duplication \\
			--thread ${task.cpus} \\
			--json /dev/null --html /dev/null > $results/ALL_REPORTS/FASTQ/TRIMMING/${sample}_SE.fastp.log
	"""
}

process paired_trimming {
	//Remove adaptaters, merge reads that can be merged
	tag "$sample"
	label 'fastp'

	input:
	tuple val(sample), file(f1), file(f2)
	val(results)

	output:
	tuple val(sample), file("${sample}_1.trimmed.fastq.gz"), file("${sample}_2.trimmed.fastq.gz")

	when:
	!params.dry && !params.help && (!mf.checkFile("$results/FASTQ/TRIMMED", sample, "q.gz") && \
	!mf.checkFile("$results/BAM/RAW", sample, "bam") && \
	!mf.checkFile("$results/BAM/FILTERED", sample, "bam") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('TRIM', params.FORCE)

	script:
	"""
	fastp -i $f1 -I $f2 --out1 ${sample}_1.trimmed.fastq.gz --out2 ${sample}_2.trimmed.fastq.gz \\
			-q 25 -l 20 -e 20 --dont_eval_duplication \\
			--thread ${task.cpus} \\
			--json /dev/null --html /dev/null > $results/ALL_REPORTS/FASTQ/TRIMMING/${sample}_PE.fastp.log
	"""
}

process taxoFilterPaired {
	tag "$sample"
	label 'kraken2'

	input:
	tuple val(sample), file(f1), file(f2)
	path krakenDB_path
	val results

	output:
	tuple val(sample), file("${sample}.classified_1.fastq.gz"), file("${sample}.classified_2.fastq.gz")

	when:
	!params.dry && !params.help && (!mf.checkFile("$results/FASTQ/TRIMMED", sample, "q.gz") && \
	!mf.checkFile("$results/BAM/RAW", sample, "bam") && \
	!mf.checkFile("$results/BAM/FILTERED", sample, "bam") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('TRIM', params.FORCE)

	script:
	"""
	mv $f1 ${sample}.classified_1.fastq.gz
	mv $f1 ${sample}.classified_2.fastq.gz
	#kraken2 --db $krakenDB_path --quick --paired $f1 $f2 --use-names \
	#	--classified-out ${sample}.classified#.fastq --unclassified-out ${sample}.unclassified#.fastq \
	#	--report ${sample}.taxo.log --thread ${task.cpus}

	#more efficient compression ? rebuilt a container with pigz ?
	#gzip -1 ${sample}.classified_1.fastq; gzip -1 ${sample}.classified_2.fastq
	#gzip -1 ${sample}.unclassified_1.fastq; gzip -1 ${sample}.unclassified_2.fastq

	#moves fastq screen report to the report folder for aggregation with multiqc
	#cp ${sample}.taxo.log $results/ALL_REPORTS/FASTQ/TAXOFILTER
	"""
}

process taxoFilterSingle {
	tag "$sample"
	label 'kraken2'

	input:
	tuple val(sample), file(f1)
	path krakenDB_path
	val results

	output:
	tuple val(sample), file("${sample}.classified.fastq.gz")

	when:
	!params.dry && !params.help && (!mf.checkFile("$results/FASTQ/TRIMMED", sample, "q.gz") && \
	!mf.checkFile("$results/BAM/RAW", sample, "bam") && \
	!mf.checkFile("$results/BAM/FILTERED", sample, "bam") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('TRIM', params.FORCE)

	script:
	"""
	kraken2 --db $krakenDB_path --quick --single $f1 --use-names \
		--classified-out ${sample}.classified.fastq --unclassified-out ${sample}.unclassified.fastq \
		--report ${sample}.taxo.log --thread ${task.cpus}
	gzip -1 ${sample}.classified.fastq; gzip -1 ${sample}.unclassified.fastq

	#moves fastq screen report to the report folder for aggregation with multiqc
	cp ${sample}.taxo.log $results/ALL_REPORTS/FASTQ/TAXOFILTER
	"""
}

process pairedNormFastq {
	tag "$sample"
	label 'bbnorm'

	input:
	tuple val(sample), file(f1), file(f2)
	val results

	output:
	tuple val(sample), file("${sample}_1.filtered.fastq.gz"), file("${sample}_2.filtered.fastq.gz")

	script:
	if(params.depth_norm > 0)
	"""
	bbnorm.sh -Xmx3g \
		target=${params.depth_norm} maxdepth=${params.depth_norm} \
		in=$f1 in2=$f2 out=${sample}_1.filtered.fastq.gz out2=${sample}_2.filtered.fastq.gz

	rm -rf $results/FASTQ/TRIMMED/${sample}_1.filtered.fastq.gz $results/FASTQ/TRIMMED/${sample}_2.filtered.fastq.gz
	mv ${sample}_1.filtered.fastq.gz $results/FASTQ/TRIMMED/${sample}_1.filtered.fastq.gz ; ln -s $results/FASTQ/TRIMMED/${sample}_1.filtered.fastq.gz ${sample}_1.filtered.fastq.gz
	mv ${sample}_2.filtered.fastq.gz $results/FASTQ/TRIMMED/${sample}_2.filtered.fastq.gz ; ln -s $results/FASTQ/TRIMMED/${sample}_2.filtered.fastq.gz ${sample}_2.filtered.fastq.gz

	"""
	else
	"""
	#not the cleanest way but simple
	ln -s $f1 ${sample}_1.normed.fastq.gz
	ln -s $f2 ${sample}_2.normed.fastq.gz
	"""
}

process singleNormFastq {
	tag "$sample"
	label 'bbnorm'

	input:
	tuple val(sample), file(f1)
	val results

	output:
	tuple val(sample), file("${sample}.filtered.fastq.gz")

	script:
	if(params.depth_norm > 0)
	"""
	bbnorm.sh -Xmx3g \
		target=${params.depth_norm} maxdepth=${params.depth_norm} \
		in=$f1 out=${sample}.filtered.fastq.gz

	rm -rf $results/FASTQ/TRIMMED/${sample}*.fastq.gz
	mv ${sample}.filtered.fastq.gz $results/FASTQ/TRIMMED/${sample}.filtered.fastq.gz
	ln -s $results/FASTQ/TRIMMED/${sample}.filtered.fastq.gz ${sample}.filtered.fastq.gz
	"""
	else
	"""
	ln -s $f1 ${sample}.filtered.fastq.gz

	rm -rf $results/FASTQ/TRIMMED/${sample}*.fastq.gz
	mv ${sample}.filtered.fastq.gz $results/FASTQ/TRIMMED/${sample}.filtered.fastq.gz
	ln -s $results/FASTQ/TRIMMED/${sample}.filtered.fastq.gz ${sample}.filtered.fastq.gz
	"""
}

process competitiveMapping {
	//screen fastq by competitive mapping of a subset of reads
	//Note : genomes reads are mapped can be changed in the data/indexed_genome by: 
	//1) adding a new reference genome (fasta) 2) indexing it by bowtie2 3) Modifying data/indexed_genome/fastq_screen.conf
	tag "$sample"
	label 'fastqscreen'

	input:
	tuple val(sample), file(f1)
	path(indexed_genome_path)
	val(results)

	when:
	!params.dry && !params.help && params.contam_check == 'compet_mapping' && \
	!mf.checkFile("$results/ALL_REPORTS/FASTQ/COMPET_MAPPING", sample, "_screen.txt")

	script:
	"""
	fastq_screen --conf $indexed_genome_path/fastq_screen.conf $f1 --aligner bowtie2 --threads ${task.cpus}

	#moves fastq screen report to the report folder for aggregation with multiqc
	cp *_screen.txt $results/ALL_REPORTS/FASTQ/COMPET_MAPPING
	"""
}

mf = new myFunctions()

workflow process_fastq {
	take: single_fastq
	take: paired_fastq
	take: index

	main:
		results = file(params.results)
		data = file(params.data)
 
		file("$data/Kraken2/Kraken16Gb").mkdirs()

		file("$results/FASTQ/TRIMMED").mkdirs()
		file("$results/FASTQ/TRIMMED").mkdirs()

		file("$results/ALL_REPORTS/FASTQ/KMERFILTER").mkdirs()
		file("$results/ALL_REPORTS/FASTQ/TRIMMING").mkdirs()
		file("$results/ALL_REPORTS/FASTQ/TAXOFILTER").mkdirs()
		file("$results/ALL_REPORTS/FASTQ/ABUNDANCE").mkdirs()
		file("$results/ALL_REPORTS/FASTQ/COMPET_MAPPING").mkdirs()

		//kmer based filter for gene oriented mode
		if(params.target_region) {
			spliceReference(index)
			paired_readCleaning(paired_fastq, spliceReference.out, results)
			single_readCleaning(single_fastq, spliceReference.out, results)
			single_trimming(single_readCleaning.out, results)
			paired_trimming(paired_readCleaning.out, results)
		} else if(params.species == 'MTBC') {
			taxoFilterPaired(paired_fastq, file(data+"/Kraken2/MTBC_DB_35"), results )
			taxoFilterSingle(single_fastq, file(data+"/Kraken2/MTBC_DB_35"), results )

			single_trimming(taxoFilterSingle.out, results)
			paired_trimming(taxoFilterPaired.out, results)
		} else {
			single_trimming(single_fastq, results)
			paired_trimming(paired_fastq, results)
		}

		//if TRIM processes forced, do not import files from results folder that otherwise could cause collisions
		if( !mf.checkFORCE('TRIM', params.FORCE) ){
			Channel.fromPath(results+"/FASTQ/TRIMMED/*.{fastq,fastq}.gz", followLinks: false)
				.map{it -> [it.simpleName.split("_1|_2")[0], it]}.groupTuple().branch{
					paired: it[1].size() == 2
					single: it[1].size() == 1
				}.set{ temp }
		} else if (mf.checkFORCE('TRIM', params.FORCE)) {
			temp = Channel.empty()
		}

		pairedNormFastq(paired_trimming.out, results)
		singleNormFastq(single_trimming.out, results)

		//merge newly trimmed fastq channel with fastq alread stored in the fastq folder into a new flux  (paired)
		all_paired_trimmed = temp.paired.map(it -> [it[0], it[1][0], it[1][1]]).mix(pairedNormFastq.out)

		//merge newly trimmed fastq channel with fastq alread stored in the fastq folder into a new flux (single end)
		all_single_trimmed = temp.single.map(it -> [it[0], it[1]]).mix(singleNormFastq.out)

		//by default, screen all fastq, can be disabled with --no_contam_check flag
		competitiveMapping(all_paired_trimmed.map(it -> [it[0], it[1]]).mix(all_single_trimmed), file(params.indexed_genome_path), results)

		emit:
			all_single_trimmed
			all_paired_trimmed
			garbage = paired_trimming.out
}
