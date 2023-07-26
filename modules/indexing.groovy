#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process bwa {
	label 'bwa'

	output:
	tuple val(name), file("${name}*")

	script:
	name = params.referenceName
	referenceSequence=file(params.referenceSequence)
	"""
	bwa-mem2 index $referenceSequence -p ${name}
	"""
}

process samtools {
	errorStrategy = "finish"
	label = 'samtools'

	output:
	tuple file("${name}.fa"), file("${name}.fa.fai")

	script:
	name = params.referenceName
	referenceSequence = file(params.referenceSequence)
	referenceAnnotation = file(params.referenceAnnotation)
	"""
	if (grep -q "gz" $referenceSequence ) ; then
		gunzip -c ${referenceSequence} > ${name}.fa
	else
		ln -s $referenceSequence ${name}.fa
	fi

	samtools faidx ${name}.fa
	"""
}


process picard {
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

process snpEff {
	label 'snpeff'
	errorStrategy='finish'

	output:
	val(true)

	when:
	!params.noGFF

	script:
	name = params.referenceName
	referenceSequence = file(params.referenceSequence)
	referenceAnnotation = file(params.referenceAnnotation)
	snpeffConfig = file(params.snpeffConfig)
	"""
	mkdir -p ${projectDir}/data/snpeff/data/${name}
	mkdir -p ${projectDir}/data/snpeff/data/genomes

	if (grep -q "gz" $referenceSequence ) ; then
		gunzip -c $referenceSequence > ${projectDir}/data/snpeff/data/genomes/${name}.fa
	else
		cp $referenceSequence ${projectDir}/data/snpeff/data/genomes/${name}.fa
	fi

	if (grep -q "gz" $referenceAnnotation ) ; then
		gunzip -c $referenceAnnotation > ${projectDir}/data/snpeff/data/${name}/genes.gff
	else
		cp $referenceSequence ${projectDir}/data/snpeff/data/${name}/genes.gff
	fi

	chromosome=\$(grep ">" ${projectDir}/data/snpeff/data/genomes/${name}.fa | tr -d '>' | cut -f1 -d ' ')

	cat $snpeffConfig > ${projectDir}/data/snpeff/config
	echo -e "${name}.genome = Mycobacterium_tuberculosis" >> ${projectDir}/data/snpeff/config
	echo -e "$name.\${chromosome}.codonTable = Mycobacterium" >> ${projectDir}/data/snpeff/config

	snpEff build -c ${projectDir}/data/snpeff/config -v ${name} -nodownload -d
	"""
}

process setExecRights {
	errorStrategy = "finish"
	label 'script'

	output:
	val(true)

	script:
	"""
	chmod +750 $baseDir/bin/*
	"""
}


workflow index {

	main:
		bwa()
		samtools()
		picard(samtools.out)

		snpEff()
		setExecRights()

	emit:
		bwa = bwa.out
		samtools_picard = picard.out
		snpeff_emit_signal = snpEff.out
		binExec_emit_signal = setExecRights.out
}
