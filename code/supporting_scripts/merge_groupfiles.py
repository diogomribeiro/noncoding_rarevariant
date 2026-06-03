#/usr/bin/python3
# Merge 3 SAIGE-GENE groupfile into one

import sys

inFile1 = open(sys.argv[1],"r")
inFile2 = open(sys.argv[2],"r")
inFile3 = open(sys.argv[3],"r")
outFile = open(sys.argv[4],"w")

## Groupfile format
# ENSG00000261456 var chr10:47880:G:A chr10:47586:C:A chr10:47925:C:T chr10:49236:C:T chr10:47523:G:T chr10:47457:G:A chr10:47548:G:A
# ENSG00000261456 anno CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS CDS

geneVars = {}
for line in inFile1:
    line = line.strip()
    spl = line.split(" ")
    if spl[1] == "var":
        gene = spl[0]
        if gene not in geneVars:
            geneVars[gene] = set()
        for var in spl[2:]:
            geneVars[gene].add(var)

for line in inFile2:
    line = line.strip()
    spl = line.split(" ")
    if spl[1] == "var":
        gene = spl[0]
        if gene not in geneVars:
            geneVars[gene] = set()
        for var in spl[2:]:
            geneVars[gene].add(var)

for line in inFile3:
    line = line.strip()
    spl = line.split(" ")
    if spl[1] == "var":
        gene = spl[0]
        if gene not in geneVars:
            geneVars[gene] = set()
        for var in spl[2:]:
            geneVars[gene].add(var)


for gene in sorted(geneVars):
    text = "%s var %s\n" % (gene, " ".join(sorted(geneVars[gene])))
    outFile.write(text)
    ncs = ["NC" for i in range(len(geneVars[gene]))]
    text = "%s annot %s\n" % (gene, " ".join(ncs))
    outFile.write(text)

print("Wrote %s genes" % len(geneVars))
outFile.close()

