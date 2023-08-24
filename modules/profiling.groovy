#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process profile {
	label 'tbdetective'

	input:
	file(vcf)
	val(results)

	output:
	file("all_strains_info.txt")

	when:
	params.taxonomy == 'barcode'

	script:
	lineageSNP = file(params.lineage_db)
	antibioSNP = file(params.antibioresistance_db)
	"""
	TB-detective -i *.vcf.gz -lin $lineageSNP -ab $antibioSNP 2> /dev/null > all_strains_info.txt
	cp all_strains_info.txt $results/all_strains_info.txt
	"""
}

process vcf2fasta {
	tag "$sample"
	label 'bcftools'

	input:
	tuple val(sample), file(vcf)
	tuple file(fasta), file(fai), file(dict)

	output:
	tuple val(sample), file("${sample}.fasta")

	when:
	params.taxonomy == 'placement'

	script:
	"""
	bcftools index $vcf
	bcftools view $vcf -V indels | bcftools view -v snps,mnps -M2 -i 'FILTER="PASS"' -Oz -o ${sample}.onlySNP.vcf.gz

	bcftools index ${sample}.onlySNP.vcf.gz
	echo ">"$sample > ${sample}.fasta
	bcftools consensus ${sample}.onlySNP.vcf.gz --fasta-ref $fasta --haplotype 2 -H SA \
		--absent '-' --missing  '-' | \
		grep -v '>' | tr -d '\n' >> ${sample}.fasta
	echo "" >> ${sample}.fasta
	"""
}

process placement {
	tag "$sample"
	label 'epa'

	input:
	tuple val(sample), file(QUERY)
	file(MSA)
	file(TREE)
	file(MODEL)
	val(results)

	output:
	tuple val(sample), file("${sample}.jplace")

	script:
	"""
	epa-ng --tree ${TREE} --ref-msa $MSA --query $QUERY -T ${task.cpus} --model $MODEL
	mv epa_result.jplace ${sample}.jplace
	cp ${sample}.jplace $results/TAXO/PLACEMENT/
	"""
}

process taxoAssignement {
	tag "$sample"
	label 'script'

	input:
	tuple val(sample), file(PLACEMENT)
	file(TAXO)

	output:
	file("${sample}_taxo_assignation.txt")

	script:
	"""
	gappa examine assign --jplace-path . --file-prefix ${sample}_placement. \
		--threads ${task.cpus} --taxon-file $TAXO --best-hit

	echo -ne "$sample\t" > ${sample}_taxo_assignation.txt
	tail -n1 ${sample}_placement.profile.tsv >> ${sample}_taxo_assignation.txt
	"""
}

process resumeTaxo {
	label 'script'

	input:
	file(assignement)
	file(results)

	script:
	"""
	echo -e "sample\tLWR\tfract\taLWR\tafract\ttaxopath" > all_phylo_placement.tsv
	cat $assignement >> all_phylo_placement.tsv
	cp all_phylo_placement.tsv $results/all_phylo_placement.tsv
	"""
}

mf = new myFunctions()

workflow profiling {
	take: annotated_vcf
	take: binExec_emit_signal
	take: index

	main:
		results = file(params.results)
		file("$results/ALL_REPORTS").mkdirs()

		file("$results/TAXO/PLACEMENT").mkdirs()
		file("$results/TAXO/LINEAGE").mkdirs()
		file("$results/TAXO/SNP_LINEAGE").mkdirs()
		file("$results/AMR").mkdirs()
		file("$results/TAXO/LINEAGE_FULL_BARCODE").mkdirs()

		strain_info = Channel.empty()
		if(params.species == 'MTBC') {
			println "Samples are assumed to belong to the M. tuberculosis Complex !"

			profile(annotated_vcf.map{it -> it[1]}.collect(), results)
			strain_info = profile.out.ifEmpty(true)
		}

		fasta = Channel.empty()
		if(params.species == 'MTBC' || params.taxonomy == 'placement') {
			vcf2fasta(annotated_vcf, index)
			placement(vcf2fasta.out, file(params.referenceMSA), file(params.referenceTree), file(params.referenceModel), file(results))
			taxoAssignement(placement.out, file(params.referenceTaxo))
			resumeTaxo(taxoAssignement.out.collect(), results)

			fasta = vcf2fasta.out
		}

	emit:
		strain_info
		fasta
}
