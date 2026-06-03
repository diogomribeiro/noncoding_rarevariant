#/usr/bin/python3
# Convert SAIGE-GENE groupfile into REGENIE annotation file for gene-based tests

import sys

inFile = open(sys.argv[1],"r")
outFile = open(sys.argv[2]+".annotations","w")
outFile2 = open(sys.argv[2]+".setlist","w")

## Input format
# ENSG00000261456 var chr10:47880:G:A chr10:47586:C:A chr10:47925:C:T chr10:49236:C:T chr10:47523:G:T chr10:47457:G:A chr10:47548:G:A
# ENSG00000261456 anno CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS

## output format 1
# 1:55039839:T:C PCSK9 LoF
# 1:55039842:G:A PCSK9 missense

## output format 2
# A1BG 19  58346922  19:58346922:C:A,19:58346924:G:A,...

geneVars = {}
varAnnot = {}
for line in inFile:
	line = line.strip()
	spl = line.split(" ")
	gene = spl[0]
	if gene not in geneVars:
		geneVars[gene] = []
		varAnnot[gene] = []
	if spl[1] == "var":
		for var in spl[2:]:
			geneVars[gene].append(var)

		chro,pos,_,_ = spl[2].split(":")
		text = "%s\t%s\t%s\t%s\n" % (gene,chro,pos,",".join(spl[2:]))
		outFile2.write(text)
	else:
		for ann in spl[2:]:
			varAnnot[gene].append(ann)

print("Genes read:", len(geneVars))

count = 0
for gene in geneVars:
	for i in range(len(geneVars[gene])):
		assert(len(geneVars[gene]) == len(varAnnot[gene]))
		count+=1
		text = "%s\t%s\t%s\n" % (geneVars[gene][i],gene,varAnnot[gene][i])
		outFile.write(text)

print("Wrote %s variants" % count)
outFile.close()
outFile2.close()