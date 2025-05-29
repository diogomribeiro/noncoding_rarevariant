#02-Feb-2023 Diogo Ribeio @ UNIL
# Script to create groupfiles for SAIGE-GENE+ input

import sys
import pandas
import numpy as np

bedFile = open(sys.argv[1],"r")
sitesFile = sys.argv[2]
outFile = open(sys.argv[3],"w")

if len(sys.argv) > 4:
    # Example format, TSV, CADD
    # chr1 12356 A C
    filterData = pandas.read_csv(sys.argv[4], sep = "\t", header = None)
    print("Read filter file 1 with %s variants" % len(filterData))
else:
    filterData = pandas.DataFrame()

### Read sites file
## Example input format, no header, TSV
# chr20	215151	C	A
sitesData = pandas.read_csv(sitesFile, sep = "\t", header = None)
sitesData = sitesData.astype({1:'int'})
print("Read %s input variants" % len(sitesData))

## If only 1 file provided, use only that file, if 2 files provided get the union of them
if len(filterData) > 0:
    sitesData = pandas.merge(sitesData, filterData, how ='inner', on = [0,1,2,3])
    sitesData.drop_duplicates()
    sitesData.reset_index(drop=True, inplace=True)

print("After filter: %s input variants" % len(sitesData))

### Read bed file
## Example input format, no header, TSV
# chr20 12315 14567 ENSG000313556 exon
totCount = 0
count = 0
geneVars = {}
for line in bedFile:
    count+=1
    if count % 1000 == 0: print("read %s bed lines.." % count)

    line = line.strip()
    spl = line.split("\t")
    chro = spl[0]
    start = int(spl[1])
    end = int(spl[2])
    gene = spl[3]
    group = spl[4]

    # Efficient implementation
    wantedVars = sitesData.loc[np.where((sitesData[0].values == chro) & (sitesData[1].values >= start) & (sitesData[1].values <= end))]

    if len(wantedVars) > 0:
        if gene not in geneVars:
            geneVars[gene] = {}
        if group not in geneVars[gene]:
            geneVars[gene][group] = set()

        for _,var in wantedVars.iterrows():
            txt = "%s:%s:%s:%s" % (var[0],var[1],var[2],var[3])
            geneVars[gene][group].add(txt)
            totCount+=1

print("%s genes with data" % len(geneVars))
print("%s total variants mapped to genes" % totCount)


## write to file
## Wanted output format (space-separated)
#ENSG00000141956 var chr21:41801250:T:C chr21:41801251:A:G
#ENSG00000141956 anno lof missense
outCount = 0
for gene in geneVars:
    text1 = "%s var" % (gene)
    text2 = "%s anno" % (gene)
    for group in geneVars[gene]:
        for var in geneVars[gene][group]:
            text1+=" %s" % var
            text2+=" %s" % group
            outCount+=1

    outFile.write(text1 + "\n")
    outFile.write(text2 + "\n")


outFile.close()
print("Wrote %s variants" % (outCount))
