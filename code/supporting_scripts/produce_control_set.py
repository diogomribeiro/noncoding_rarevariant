#!/usr/bin/python3
import sys
import pandas as pd

window = 1000000
seed = 999

## INPUT
numberVarsFile = open(sys.argv[1],"r") # gene\tnumber
listPosExcludeFile = sys.argv[2] # all_vars_JAVIERRE_ABC_CRD_CDS.txt.gz --> this is only CHROM and POS
listVarsFile = open(sys.argv[3],"r") # groupfile format, list of variants passing CADD 15 in cis-window of gene
outFile = open(sys.argv[4],"w") # groupfile format

numberVars = {} # key -> gene, val -> number of vars to subsample
for line in numberVarsFile:
	spl = line.strip().split("\t")
	gene = spl[0]
	variants = int(spl[1])
	numberVars[gene] = variants

print("read number vars for %s genes" % len(numberVars))

listPosExclude = pd.read_csv(listPosExcludeFile, sep = "\t", header = None, names=["chrom","pos"], dtype={'chrom' : str, 'pos' : int})

print("read %s positions to exclude" % len(listPosExclude))

count = 0
totCount = 0
for line in listVarsFile:
	totCount+=1
	line = line.strip()
	spl = line.split(" ")
	if spl[1] != "var":
		continue
	gene = spl[0]
	variants = spl[2:]
	if gene not in numberVars:
		# Meaning: no non-coding variants were present, thus no control needed for this gene
		continue

	### Exclude set of variants
	data = [item.split(':') for item in variants]
	columns = ['chrom', 'pos', 'a1', 'a2']
	variantsDF = pd.DataFrame(data, columns=columns)
	variantsDF['pos'] = variantsDF['pos'].astype(int)

	# speed up search to gene window
	genePos = int(variants[0].split(":")[1])
	genePosMin = genePos - window * 2
	genePosMax = genePos + window * 2
	possibleExclude = listPosExclude[listPosExclude["pos"] > genePosMin][listPosExclude["pos"] < genePosMax]

	mergedData = pd.merge(variantsDF, possibleExclude, on=['chrom', 'pos'], how='outer', indicator=True)
	varPool = mergedData[mergedData['_merge'] == 'left_only']
	varPool = varPool.drop(columns=['_merge']).reset_index(drop=True)

	assert(varPool.shape[0] <= variantsDF.shape[0])

	### Subsample variants to match number from non-coding analysis
	nVar = numberVars[gene]
	if nVar > varPool.shape[0]:
		print("PROBLEM for gene %s: not enough variants to subsample" % gene)
		continue

	subsampledDF = varPool.sample(n=nVar, random_state = seed, axis = 0)
	subsampledDF['tag'] = subsampledDF.apply(lambda row: ':'.join(map(str, row)), axis=1)
	finalVars = {row for row in subsampledDF['tag'].values}

	assert(subsampledDF.shape[0] == nVar)

	count+=1
	if count % 1000 == 0:
		print("wrote %s genes.." % count)

	### write groupfile as output 
	outFile.write("%s var %s\n" % (gene," ".join(sorted(finalVars))))
	outFile.write("%s anno %s\n" % (gene," ".join(["CONTROL"] * len(finalVars))))

print("Wrote variants for %s genes, out of %s in input" % (count, totCount))

outFile.close()
print("FINISHED!")
