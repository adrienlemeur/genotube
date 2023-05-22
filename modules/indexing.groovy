#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process bwa_index {
	label 'bwa'

	input:
	file(referenceSequence)

	output:
	tuple val(referenceSequence.simpleName), file("${referenceSequence.simpleName}*")

	when:
	!params.help && !params.dry

	script:
	"""
	bwa-mem2 index $referenceSequence -p $referenceSequence.simpleName
	"""
}

process samtools_index {
	errorStrategy = "finish"
	label = 'samtools'

	input:
	file(fasta)

	output:
	tuple file("${fasta.simpleName}.fa"), file("${fasta.simpleName}.fa.fai")

	when:
	!params.help && !params.dry

	script:
	"""
	{ 
		zcat $fasta > ${fasta.simpleName}.fa
	} || {
		mv $fasta ${fasta.simpleName}.fa
       }
	samtools faidx ${fasta.simpleName}.fa
	"""
}


process picard_index {
	errorStrategy = "finish"
	label 'GATK'

	input:
	tuple file(fasta), file(fai)

	output:
	tuple file(fasta), file(fai), file("*.dict")

	script:
	"""
	gatk CreateSequenceDictionary -R $fasta
	"""
}

process buildSnpeffDB {
	label 'snpeff'
	errorStrategy='finish'

	input:
	file(snpeffConfig_file)
	file(referenceSequence)
	file(referenceAnnotation)

	output:
	val(true)

	when:
	!params.help && !params.dry && !params.noGFF

	script:
	"""
	path_to_config=\$(dirname \$(realpath snpEff.config))

	mkdir -p \$path_to_config/data/${referenceSequence.simpleName}
	mkdir -p \$path_to_config/data/genomes

	zcat $referenceSequence > \$path_to_config/data/genomes/${referenceSequence.simpleName}.fa
	zcat $referenceAnnotation > \$path_to_config/data/${referenceSequence.simpleName}/genes.gff

	echo "tteeeeeeeeest"
	echo "$referenceSequence"
	echo "$referenceSequence.simpleName"
	echo "$referenceAnnotation"
	
	snpEff build -c $snpeffConfig_file -v ${referenceSequence.simpleName} -nodownload -d
	"""
}

process binExecutionRights {
	errorStrategy = "finish"
	label 'script'

	input:
	file(binPath)

	output:
	val(true)

	when:
	!params.help

	script:
	"""
	chmod +777 $binPath/*
	"""
}


workflow index {

	main:
		snpeffConfig = file("data/snpeff/snpEff.config")

		bwa_index(file(params.referenceSequence))
		samtools_index(file(params.referenceSequence))
		picard_index(samtools_index.out)

		buildSnpeffDB(snpeffConfig, file(params.referenceSequence), file(params.referenceAnnotation))
		binExecutionRights(file("bin"))

	emit:
		bwa_index = bwa_index.out
		samtools_picard_index = picard_index.out
		snpeff_emit_signal = buildSnpeffDB.out
		binExec_emit_signal = binExecutionRights.out
}
