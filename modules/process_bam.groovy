#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process qualityCheck {
	tag "$sample"
	label 'bedtools'

	input:
	tuple val(sample), file(bam)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), env(covdepth), env(cov), file(bam), emit: file
	tuple val(sample), env(covdepth), env(cov), emit: covinfo
	val true, emit: signal

	when:
	(!mf.checkFile("$results/BAM/FILTERED", sample, "bam") && \
	!mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('MAP', params.FORCE)

	script:
	if(params.target_region){region = file(params.target_region)} else {region = false}
	if( params.target_region == false )
		"""
		covdepth=\$(bedtools genomecov -ibam $bam | grep genome | LC_ALL=C awk 'NR>1 {sum+=\$2*\$3} END {print sum/\$4}' | cut -f1 -d '.')
		cov=\$(bedtools genomecov -ibam $bam -max 1 | LC_ALL=C awk 'NR==2 {print \$5*100}' | cut -d '.' -f1)

		echo -e "$sample\t\$covdepth\t\$cov" > ${sample}_QC.txt
		cat ${sample}_QC.txt > $results/ALL_REPORTS/BAM/QUAL/${sample}_QC.txt

		"""
	else if( params.target_region != false ){
		"""
		covdepth=\$(bedtools coverage -a $region -b $bam -mean | LC_ALL=C awk '{total+=\$4} END {print total/NR}' | cut -f1 -d '.')
		cov=\$(bedtools coverage -a $region -b $bam | LC_ALL=C awk '{total+=\$6 ; base_covered+=(\$6*\$7)} END {print 100*(base_covered/total)}' | cut -f1 -d '.')

		echo -e "$sample\t\$covdepth\t\$cov" > ${sample}_QC.txt
		cat ${sample}_QC.txt > $results/ALL_REPORTS/BAM/QUAL/${sample}_QC.txt
		"""
	} else
		error "You should not be there"
}

process extQualityCheck {
	tag "$sample"
	label 'bedtools'

	input:
	tuple val(sample), file(bam)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), env(covdepth), env(cov), file(bam), emit: file
	tuple val(sample), env(covdepth), env(cov), emit: covinfo
	val true, emit: signal

	when:
	(!mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") && \
	!mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('CALL', params.FORCE)

	script:
	if(params.target_region){region = file(params.target_region)} else {region = false}
	if( params.target_region == false )
		"""
		covdepth=`bedtools genomecov -ibam $bam | grep genome | LC_ALL=C awk 'NR>1 {sum+=\$2*\$3} END {print sum/\$4}' | cut -f1 -d '.'`
		cov=`bedtools genomecov -ibam $bam -max 1 | LC_ALL=C awk 'NR==2 {print \$5*100}' | cut -d '.' -f1`

		echo -e "$sample\t\$covdepth\t\$cov" > ${sample}_QC.txt
		cat ${sample}_QC.txt > $results/ALL_REPORTS/BAM/QUAL/${sample}_QC.txt
		"""
	else if( params.target_region != false ){
		region = file(params.target_region)
		"""
		covdepth=`bedtools coverage -a $region -b $bam -mean | LC_ALL=C awk '{total+=\$4} END {print total/NR}' | cut -f1 -d '.'`
		cov=`bedtools coverage -a $region -b $bam | LC_ALL=C awk '{total+=\$6 ; base_covered+=(\$6*\$7)} END {print 100*(base_covered/total)}' | cut -f1 -d '.'`

		echo -e "$sample\t\$covdepth\t\$cov" > ${sample}_QC.txt
		cat ${sample}_QC.txt > $results/ALL_REPORTS/BAM/QUAL/${sample}_QC.txt
		"""
	} else
		error "You should not be there"
}

process resumeCoverage {

	input:
	val(all_files)
	val(results)

	script:
	"""
	echo -e "sample\tsequencing_depth\tcoverage" > $results/samples_coverage.txt
	cat \$(find $results -name "*_QC.txt") >> $results/samples_coverage.txt
	"""
}

process downSamplingAlignedRead {
	tag "$sample"
	label "samtools"

	input:
	tuple val(sample), file(bam)
	val(results)

	output:
	tuple val(sample), file("${sample}.subsample.fastq.gz")

	when:
	!mf.checkFile("$results/ALL_REPORTS/BAM/BRACKEN", sample, ".abundance.log")

	script:
	"""
	frac=\$( samtools idxstats $bam | cut -f3 | awk 'BEGIN {total=0} {total += \$1} END {frac=200000/total; if (frac > 1) {print "1.00"} else {print frac}}' )
	samtools view -bs \$frac $bam | samtools fastq > ${sample}.subsample.fastq
	bgzip ${sample}.subsample.fastq
	"""
}

process downloadKrakenDB {
	label 'script'

	input:
	val(data)

	output:
	val true

	when:
	( (params.prebuilt_K2_DB == '8G' || params.prebuilt_K2_DB == '16G') && \
	!file("$data/Kraken2/${params.prebuilt_K2_DB}/inspect.txt").exists() ) || custom_K2_DB

	script:
	println "Downloading the Kraken2 database, might take a while depending on your connection (< ${params.prebuilt_K2_DB})"
	if( params.prebuilt_K2_DB == '8G' ){
		database="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20221209.tar.gz"
	} else if ( params.prebuilt_K2_DB == '16G' ) {
		database="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20220926.tar.gz"
	}
	db = params.prebuilt_K2_DB
	"""
	mkdir $data/Kraken2/${db}
	wget ${database} -t ${task.cpus} -O ${db}.tar.gz

	gunzip ${db}.tar.gz
	tar -xf ${db}.tar --directory $data/Kraken2/${db}
	touch $data/Kraken2/${db}/inspect.txt
	"""
}

process taxoClass {
	tag "$sample"
	label 'kraken2'

	input:
	val downloaded
	tuple val(sample), file(f1)
	val data

	output:
	tuple val(sample), file("${sample}.taxo.log")

	script:
	if(params.custom_K2_DB){ DB = file(params.custom_K2_DB) }
	if(params.custom_K2_DB)
	"""
	kraken2 --db $DB --single $f1 --use-names \
		--report ${sample}.taxo.log --thread ${task.cpus}
	"""
	else
	"""
	kraken2 --db $data/Kraken2/${params.prebuilt_K2_DB}/ --single $f1 --use-names \
		--report ${sample}.taxo.log --thread ${task.cpus}
	"""
}
	
process abundance {
	tag "$sample"
	label 'bracken'

	input:
	tuple val(sample), file(kraken_log)
	val krakenDB_path
	val results

	script:
	"""
	#lower detection treshold
	bracken -d $krakenDB_path -i $kraken_log -o ${sample}.bracken.log -w ${sample}.abundance.log
	cp ${sample}.abundance.log $results
	"""
}

process deduplication {
	tag "$sample"
	label 'GATK_SPARK'

	input:
	tuple val(sample), file(bam)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}.dedup.bam"), file("*.bai")

	when:
	((! mf.checkFile("$results/BAM/FILTERED", sample, "bam") && \
	! mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") && \
	! mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('CALL', params.FORCE) ) && !params.target_region

	script:
	"""
	gatk MarkDuplicatesSpark \\
		-I $bam --create-output-bam-index \\
		-M $results/ALL_REPORTS/BAM/DEDUP/${sample}_dedup_metrics.txt \\
		-O ${sample}.dedup.bam \\
		--remove-all-duplicates

	rm -rf $results/BAM/FILTERED/${sample}.dedup.bam

	mv ${sample}.dedup.bam $results/BAM/FILTERED/${sample}.dedup.bam
	ln -s $results/BAM/FILTERED/${sample}.dedup.bam ${sample}.dedup.bam
	"""
}

process statFastQC {
	tag "$sample"
	label 'fastqc'

	input:
	tuple val(sample), file(fastq)
	val(results)

	output:
	val true

	when:
	false

	script:
	"""
	fastqc -t ${task.cpus} --extract $fastq -o ${results}/ALL_REPORTS/BAM/QUAL/
	"""
}

process indexBam {
	tag "$sample"
	label 'samtools'

	input:
	tuple val(sample), file(bam)

	output:
	tuple val(sample), file(bam), file("*.bai")

	when:
	((! mf.checkFile("$results/BAM/FILTERED", sample, "bam") && \
	! mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") && \
	! mf.checkFile("$results/VCF/RAW", sample, "vcf.gz") ) || mf.checkFORCE('CALL', params.FORCE) )

	script:
	"""
	samtools index $bam
	"""
}

mf = new myFunctions()

workflow process_bam {
	take: all_bam
	take: index

	main:
		results = file(params.results)
		data = file(params.data)

		file("$results/ALL_REPORTS/BAM").mkdirs()
		file("$results/ALL_REPORTS/BAM/DEDUP").mkdirs()
		file("$results/ALL_REPORTS/BAM/QUAL").mkdirs()

		file("$results/BAM/RAW").mkdirs()
		file("$results/BAM/FILTERED").mkdirs()

		file("$results/VCF/RAW").mkdirs()
   		file("$results/DELETION_REGION").mkdirs()

		//Quality check
		if(! mf.checkFORCE('COVERAGE_CHECK', params.SKIP) ) {
			qualityCheck(all_bam, index, results)
			//Discard sample with sequencing depth <= 10 && coverage >= 80%
			qualityCheck.out.file.filter{it[1].toInteger() >= 10 && it[2].toInteger() >= 80}.map{it -> [ it[0], it[3] ] }.set{good_bam}
		} else {
			//no quality check, pass all bam downstream
			all_bam.set{ good_bam }
		}

//		deduplication(good_bam, index, results)

		//deduplication is slow, add already processed bam
		if( ! mf.checkFORCE('MAP', params.FORCE) ) {
			old_processed_bam = Channel.fromPath(results+"/BAM/FILTERED/*.bam").map{it -> [it.simpleName, it]}
		} else {
			old_processed_bam = Channel.empty()
		}

		if(! mf.checkFORCE('COVERAGE_CHECK', params.SKIP) ) {	//checkForce name is confusing, should either be checkSKIP or a generic function for both since params.[FORCE/SKIP] is specified
			extQualityCheck(old_processed_bam, index, results)
			extQualityCheck.out.file.filter{it[1].toInteger() >= 10 && it[2].toInteger() >= 80}.map{it -> [ it[0], it[3] ] }.set{ good_old_bam }

			coverage_info = qualityCheck.out.covinfo.mix(extQualityCheck.out.covinfo)
			bamChecked = extQualityCheck.out.signal.collect().ifEmpty(true).mix(qualityCheck.out.signal.ifEmpty(true)).collect()
		} else {
			old_processed_bam.set{ good_old_bam }

			//when quality check is disabled, set samples coverage & depth to null
			coverage_info = good_old_bam.mix(good_bam).map{it -> [it[0]]}.flatten().unique().map{it -> [ it, null, null ]}
			bamChecked = coverage_info.collect().ifEmpty(true)
		}

		//index extern files
		indexBam(good_bam.mix(good_old_bam))

		//add them to the flow
		all_processed_bam = indexBam.out

		file("$results/ALL_REPORTS/BAM/KRAKEN").mkdirs()
		file("$results/ALL_REPORTS/BAM/BRACKEN").mkdirs()

		downloadKrakenDB(data)

		downSamplingAlignedRead(all_bam, results)
		statFastQC(downSamplingAlignedRead.out, results)
		taxoClass(downloadKrakenDB.out.ifEmpty(true), downSamplingAlignedRead.out, data)
		abundance(taxoClass.out, file(data+"/Kraken2/$params.prebuilt_K2_DB"), results+"/ALL_REPORTS/BAM/BRACKEN")

		//end of quality check triggers coverage file concatenation / should be fixed
		resumeCoverage(bamChecked.collect().ifEmpty(true), results)

	emit:
		all_processed_bam
		coverage_info
}
