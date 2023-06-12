#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process getAntibioRType {
	tag "$sample"
	label 'R'

	input:
	tuple val(sample), file(vcf)
	file(antibioResistanceList)
	val(start)
	val(results)

	output:
	val(true)

	when:
	!mf.checkFile("$results/AMR", sample, "_antibio_info.txt") || mf.checkFORCE('ANN', params.FORCE)

	script:
	"""
	zcat $vcf | grep -v FAIL | grep -v '#' > ${sample}.vcf
	isStrainResistant.sh ${sample}.vcf $antibioResistanceList > ${sample}

	getDT.R ${sample}
	mv ${sample}.refined.txt ${sample}_antibio_info.txt
	cp ${sample}_antibio_info.txt $results/AMR

	rm -rf ${sample}.vcf ${sample}
	"""
}


process getLineage {
	tag "$sample"
	label 'python'

	input:
	tuple val(sample), file(vcf)
	file(lineageSNP)
	val(results)

	output:
	val true

	when:
	!mf.checkFile("$results/TAXO/LINEAGE_FULL_BARCODE", sample, "_full_lineage_report.txt") && \
	!mf.checkFile("$results/TAXO/LINEAGE_FULL_BARCODE", sample, "_full_lineage_report.txt") && \
	params.taxonomy == 'barcode' || mf.checkFORCE('ANN', params.FORCE)

	script:
	"""
	{
		cut -f3 $lineageSNP > pos
		zcat $vcf | grep -v FAIL | grep -Fwf pos | cut -f1-5 > ${sample}.tmp
		cut -f2 ${sample}.tmp > backPos
		grep -Fwf backPos $lineageSNP > refLineage.bed

		getLineage.py -i ${sample}.tmp -R refLineage.bed > $results/TAXO/SNP_LINEAGE/${sample}_lineage_info.txt
		mv *_full_lineage_report.txt $results/TAXO/LINEAGE_FULL_BARCODE/${sample}_full_lineage_report.txt
	} || {
		touch probably_not_tuberculosis
	}
	"""
}

process resumeStrainInfo {
	label 'script'

	input:
	tuple val(sample), file(vcf)
	file(antibioProfile)
	file(lineageProfile)
	val(results)

	output:
	tuple val(sample), file(vcf), file("all_strain_info.txt")
	
	when:
	params.taxonomy == 'barcode'

	script:
	"""
	echo -ne "sample\tsignal\tlineage\tsublineage\tDR_type\t" > all_strain_info.txt
	echo -ne "amikacin\taminoglycosides\tbedaquiline\tcapreomycin\tciprofloxacin\tclofazimine\t" >> all_strain_info.txt
	echo -ne "cycloserine\tdelamanid\tethambutol\tethionamide\t" >> all_strain_info.txt
	echo -ne "fluoroquinolones\tisoniazid\tkanamycin\tlevofloxacin\tlinezolid\t" >> all_strain_info.txt
	echo -ne "moxifloxacin\tofloxacin\tpara-aminosalicylic_acid" >> all_strain_info.txt
	echo -ne "\tpyrazinamide\trifampicin\tstreptomycin\n" >> all_strain_info.txt

	join <(sort -k1 $results/TAXO/SNP_LINEAGE/*) <(sort -k1 $results/AMR/*) -e 'NA' -t \$'\t' | sort -k4 >> all_strain_info.txt
	cp all_strain_info.txt $results/all_strain_info.txt
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
	take: raw_vcf
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

			antibioResistanceList = file("data/antibioResistance.txt")

			getLineage(annotated_vcf, file(params.lineageSNP), results)
			getAntibioRType(annotated_vcf, antibioResistanceList, binExec_emit_signal, results)

			resumeStrainInfo(annotated_vcf.collect(), getAntibioRType.out.collect().ifEmpty(true), getLineage.out.collect().ifEmpty(true), results)
			strain_info = resumeStrainInfo.out.ifEmpty(true)
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
