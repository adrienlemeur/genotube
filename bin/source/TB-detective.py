#!/usr/bin/env python3

from cyvcf2 import VCF
import numpy as np
import argparse
import re, sys, gc

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(
			prog="TB-barcode",
			description="""TB-detective is a light-weight binary to identify the lineage, sublineage and antibiotic resistance of a Mycobacterium tuberculosis sample from a VCF annotated with snpEff. It was written in python with cyvcf2 and compiled with Nuitka3.""",
			epilog="Written by Adrien Le Meur, v.0.3")

parser.add_argument('-i', type=str, nargs='+', required=True, help='a single sample VCF aligned on H37Rv genome')
parser.add_argument('-lin', type=str, nargs=1, required=False, help='a tab separated table with 1-based SNP position, the ALT nucleotide and the associated lineage')
parser.add_argument('-ab', type=str, nargs=1, required=False, help='a tab separated table with genes, mutations (snpEff prot. or nuc. mutation annotation) and associated antibiotic resistance')

args = parser.parse_args()

lineage = args.lin
input_vcf = args.i
antibio = args.ab

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

#antibioresistance & lineage information are stored in two dictionaries for fast accession
amr_dict = {}
sample_antibio_dict = {}

lineage_dict = {}
for line in open(lineage[0]):
	line = line.rstrip().split("\t")
	lineage_dict[line[0]] = {}
	lineage_dict[line[0]]["VAR"] = line[1]
	lineage_dict[line[0]]["LIN"] = line[2]

for line in open(antibio[0]):
	line = line.rstrip().split("\t")
	key = "_".join( [ line[0], line[1] ] )
	amr_dict[key] = {}
	amr_dict[key]["gene"] = line[0]
	amr_dict[key]["position"] = line[1]
	amr_dict[key]["resistance"] = line[2]
	sample_antibio_dict[line[2]] = "FALSE" #additional dict with antibiotic resistance class as key and TRUE (resistant) or FALSE (sensitive)

#print the header
print("sample", "signal_type", "lineages", "sublineages", "drug_res_type", "resistance_count", "\t".join(sample_antibio_dict.keys()), sep = "\t")

#for every VCF input in -i
for sample in input_vcf:
	vcf = VCF(sample)

	sample = vcf.samples[0]
	eprint("Doing "+sample+" ...")

	sample_lineage_dict = {}

	#printing all lineage SNP found
	LBC = open("all_samples_full_barcode.txt", "w")

	antibio_count = 0

	#printing all antibiotic SNP found
	AB = open("all_samples_AB.txt", "w")

	for i in sample_antibio_dict:
		sample_antibio_dict[i] = "FALSE"

	for v in vcf:
		if v.INFO.get("ANN"):

			if(v.FILTER is None):
				FILTER = 'PASS'
			else:
				FILTER = v.FILTERS[0]

			ANN = v.INFO["ANN"].split("|")
			if(ANN[7] == "protein_coding"):
				nucleotide_key = "_".join([ANN[3],ANN[9]])
				protein_key = "_".join([ANN[3],ANN[10]])

				if(nucleotide_key in amr_dict):
					AB.write(sample+"\t"+amr_dict[nucleotide_key]["gene"]+"\t"+amr_dict[nucleotide_key]["position"]+"\t"+FILTER+"\t"+amr_dict[nucleotide_key]["resistance"]+'\n')
					if FILTER == 'PASS':
						antibio_count+=1
						sample_antibio_dict[amr_dict[nucleotide_key]["resistance"]] = "TRUE"
				if(protein_key in amr_dict):
					AB.write(sample+"\t"+amr_dict[protein_key]["gene"]+"\t"+amr_dict[protein_key]["position"]+"\t"+FILTER+"\t"+amr_dict[protein_key]["resistance"]+'\n')
					if FILTER == 'PASS':
						antibio_count+=1
						sample_antibio_dict[amr_dict[protein_key]["resistance"]] = "TRUE"
		if v.is_snp & (FILTER == "PASS"):
			variant_pos = str(v.POS)
			if(variant_pos in lineage_dict):
				if lineage_dict[variant_pos]["LIN"] not in sample_lineage_dict:
					sample_lineage_dict[lineage_dict[variant_pos]["LIN"]] = 1
				else:
					sample_lineage_dict[lineage_dict[variant_pos]["LIN"]] = sample_lineage_dict[lineage_dict[variant_pos]["LIN"]] + 1
				LBC.write(sample+"\t"+variant_pos+"\t"+v.REF+"\t"+lineage_dict[variant_pos]["VAR"]+"\t"+lineage_dict[variant_pos]["LIN"]+"\t"+FILTER+"\n")
	vcf.close()
	AB.close()
	LBC.close()

	if len(sample_lineage_dict) == 0:
		signal = "NO_SNP"
		main_sublineages = ["NONE"]
		main_lineages = ["NONE"]
	else:
		sample_lineage_list = list(sample_lineage_dict)
		sample_lineage_list.sort(reverse=True, key=lambda x: (len(x), x))
		main_sublineages = [ sample_lineage_list.pop(0) ]

		while len(sample_lineage_list) > 0:
			new = True
			first_key = sample_lineage_list[0]
			for i in main_sublineages:
				if re.search(first_key, i) is not None:
					new = False
			if new == True:
				main_sublineages.append(first_key)
			sample_lineage_list.pop(0)

			if len(main_sublineages) == 1:
				signal = "CLEAR"
			elif len(main_sublineages) == 2:
				signal = "COINFECTION"
			else:
				signal = "CONTAMINATION"
		
		main_lineages = set([ i.split('.')[0] for i in main_sublineages ])

	DR_type = "SENSITIVE"

	if("isoniazid" in sample_antibio_dict and "rifampicin" in sample_antibio_dict):
		if(sample_antibio_dict["isoniazid"] == "TRUE" and sample_antibio_dict["rifampicin"] == "TRUE"):
			DR_type = "MDR"
			if("fluoroquinolones" in sample_antibio_dict):
				if(sample_antibio_dict["fluoroquinolones"] == "TRUE"):
					if("capreomycin" in sample_antibio_dict):
						if(sample_antibio_dict["capreomycin"] == "TRUE"):
							DR_type = "XDR"
					if("kanamycin" in sample_antibio_dict):
						if(sample_antibio_dict["kanamycin"] == "TRUE"):
							DR_type = "XDR"
					if("amikacin" in sample_antibio_dict):
						if(sample_antibio_dict["amikacin"] == "TRUE"):
							DR_type = "XDR"
		elif "rifampicin" in sample_antibio_dict:
			if sample_antibio_dict["rifampicin"] == "TRUE":
				DR_type = "RR"

	print(sample, signal, ";".join(list(main_lineages)), ";".join(list(main_sublineages)), DR_type, str(antibio_count), "\t".join(sample_antibio_dict.values()), sep = "\t")
