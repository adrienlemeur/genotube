#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process vcf2Fasta {
	label 'bcftools'

	input:
	tuple val(sample), file(vcf)
	tuple file(fasta), file(fai), file(dict)
	file(excludeRegions)

	output:
	file("${sample}_rebuilt.fasta")

	when:
	params.tree_build != 'none'

	script:
	"""
	bcftools index $vcf
	bcftools filter $vcf -i "MAX(INFO/AF)>0.80 && MAX(FORMAT/AD)>=5 && QUAL>30" --SnpGap 2 --IndelGap 9 | \
		bcftools view -v snps,mnps -M2 -Oz -o ${sample}.onlySNP.vcf.gz

	bcftools index ${sample}.onlySNP.vcf.gz
	echo ">"$sample > ${sample}_rebuilt.fasta

	bcftools consensus ${sample}.onlySNP.vcf.gz --fasta-ref $fasta --haplotype 1 -H SR | \
		grep -v '>' | tr -d '\n' >> ${sample}_rebuilt.fasta
	echo -ne "\n" >> ${sample}_rebuilt.fasta
	"""
}

process treeBuild {
	label 'raxml'

	input:
	file(fasta)
	val(results)

	output:
	file("full_tree.raxml.bestTree")

	when:
	params.tree_build == 'classic'

	script:
	"""
	cat $fasta > MSA.fasta
	cat MSA.fasta > $results/MSA.fasta

	raxml-ng --all --msa MSA.fasta --model DNA --bs-trees 100 --prefix full_tree \\
		--threads ${task.cpus}
	cp -f full_tree* $results
	"""
}

process fast_treeBuild {
	label 'fasttree'

	input:
	file(fasta)
	val(results)

	output:
	tuple file("MSA.fasta"), file("full_tree.fasttree.bestTree")

	when:
	params.tree_build == 'fastest'

	script:
	"""
	cat $fasta > MSA.fasta
	cat MSA.fasta > $results/MSA.fasta

	goalign subsites -i MSA.fasta --informative > MSA_informative_sites.fasta

	export OMP_NUM_THREADS=${task.cpus}
	FastTree -gtr -nt MSA_informative_sites.fasta > full_tree.fasttree.bestTree
	cp -f full_tree.fasttree.bestTree $results
	"""
}

process brenlenRecalibration {
	label 'raxml'

	input:
	tuple file(msa), file(tree)
	val(results)

	output:
	file("full_tree.raxml.bestTree")

	script:
	"""
	gotree resolve -i $tree > full_tree.resolved.bestTree
	raxml-ng --evaluate --msa $msa --tree full_tree.resolved.bestTree --model DNA --prefix full_tree --brlen scaled --threads ${task.cpus}
	cp -f full_tree* $results
	"""
}

process treeAnnotation {
	label 'pastml'

	input:
	file(tree)
	file(strain_info)
	val(outgroup)
	val(results)

	when:
	params.tree_build != 'none'

	script:
	if(outgroup != false){
		outgroup = file(outgroup)
		"""
		gotree reroot outgroup -i full_tree.raxml.bestTree -l $outgroup > full_tree.rerooted.raxml.bestTree
		pastml -t full_tree.rerooted.raxml.bestTree -d all_strain_info.txt -c lineage sublineage DR_type --prediction_method DOWNPASS --work_dir $results/TREE/ANN --upload_to_itol
		"""
	} else if (outgroup == false) {
		"""
		pastml -t $tree -d all_strain_info.txt -c lineage sublineage DR_type --prediction_method DOWNPASS --work_dir $results/TREE/ANN --upload_to_itol
		"""
	} else {
		error "Should not be there"
	}
}

mf = new myFunctions()

workflow build_tree {
	take: filterVariants
	take: strain_info
	take: index

	main:
		results = file(params.results)
		file("$results/TREE/RAW").mkdirs()
		file("$results/TREE/ANN").mkdirs()

		excludeRegions = file("data/problematic_H37Rv_maping_regions.bed")
		old_vcf = Channel.fromPath("$results/VCF/FILTERED/*.vcf.gz").map{it -> [it.simpleName, it]}
		allVCF = filterVariants.join(old_vcf, remainder: true).map{if(it[1] == null){[it[0], it[2] ] } else {[ it[0], it[1] ]} }

		vcf2Fasta(allVCF, index, excludeRegions)
		treeBuild(vcf2Fasta.out.collect(), results+"/TREE/RAW")
		fast_treeBuild(vcf2Fasta.out.collect(), results+"/TREE/RAW")
		brenlenRecalibration(fast_treeBuild.out,results+"/TREE/RAW")

		treeAnnotation(treeBuild.out.mix(brenlenRecalibration.out), strain_info, params.outgroup, results)
}
