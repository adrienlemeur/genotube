#!/usr/bin/env python3

import sys
vcf = sys.argv[sys.argv.index("-i")+1]
refVCF = sys.argv[sys.argv.index("-R")+1]

variants_dict = {}
with open(refVCF) as reference:
	for line in reference:
		line = line.rstrip().split("\t")
		CHROM = line[0]; POS = line[2];	LINEAGE = line[3]; VAR = line[4]

		if CHROM not in variants_dict:
			variants_dict[CHROM] = {}

		if POS not in variants_dict[CHROM]:
			variants_dict[CHROM][POS] = {}

		if VAR not in variants_dict[CHROM][POS]:
			variants_dict[CHROM][POS][VAR] = LINEAGE

sublineage_dict = {}
with open(vcf) as target:
	for line in target:
		line = line.rstrip().split("\t")
		CHROM = line[0]; POS = line[1];	VAR = line[4]

		if CHROM in variants_dict:
			if POS in variants_dict[CHROM]:
				if VAR in variants_dict[CHROM][POS]:
					if variants_dict[CHROM][POS][VAR] not in sublineage_dict:
						sublineage_dict[variants_dict[CHROM][POS][VAR]] = 1
					else:
						sublineage_dict[variants_dict[CHROM][POS][VAR]] = sublineage_dict[variants_dict[CHROM][POS][VAR]]+1

lineages = set([k.split(".")[0]  for  k in  sublineage_dict.keys()])

print(vcf.split(".")[0], end='\t')

if(len(lineages) == 1):
	print("clear", end='\t')
else:
	print("mixed", end='\t')

#SOMME DU NOMBRE DE MARQUEURS PAR LIGNEE
lineage_dict = {}
for i in sublineage_dict:
	lineage = i.split(".")[0]
	if lineage not in lineage_dict:
		lineage_dict[lineage] = {}
		lineage_dict[lineage]["count"] = 0 + sublineage_dict[i]
		lineage_dict[lineage][i] = sublineage_dict[i]
	else:
		lineage_dict[lineage]["count"] = lineage_dict[lineage]["count"] + sublineage_dict[i]
		lineage_dict[lineage][i] = sublineage_dict[i]

#MEILLEUR LIGNEE
max_lineage_marker = 0
best_lineage = []

for i in lineage_dict:
	if(lineage_dict[i]["count"] > max_lineage_marker):
		max_lineage_marker = lineage_dict[i]["count"]
		best_lineage = i


#MEILLEURE SOUS-LIGNEE
if (len(lineage_dict) > 0) :
	print(best_lineage, end = "\t")
	if best_lineage:
		print(sorted([i for i in lineage_dict[best_lineage] if i != 'count'], reverse=True)[0])
	else :
		print(i)
else :
	print("NA", "NA")

file_name = vcf.split(".")[0]+"_full_lineage_report.txt"
with open(file_name, 'a') as output:

	output.write("lineage\tsublineage\tIndependant Supporting SNP count\n")
	for i in lineages:
		for j in lineage_dict[i]:
			if(j != "count"):
				line = str(i)+"\t"+str(j)+"\t"+str(lineage_dict[i][j])+"\n"
				output.write(line)
