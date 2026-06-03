#!/usr/bin/env python3.6

import argparse
import gzip
import queue

from util.log.Logger import Logger
from util.time.Timer import Timer

#===============================================================================
DESC_COMMENT = "Script to calculate burden of rare variants in given genomic regions."
SCRIPT_NAME = "ParseCRDBlueprint.py"
#===============================================================================

"""
#===============================================================================
@author: Diogo Ribeiro
@date: 4 May 2022
@copyright: Copyright 2022, University of Lausanne
"Script to list peak regions associated with genes through CRDs.
#===============================================================================
"""

class ParseCRDBlueprint(object):
           
    def __init__(self, crdFile, geneLinkFile, outFile):

        Timer.get_instance().start_chrono()

        self.crdFile = crdFile
        self.geneLinkFile = geneLinkFile
        self.outFile = outFile


    def read_crd_file(self):
        """
        Read CRD file (format from Olivier Delaneau).

        Example format (TSV, no header):
            chr1    235531  236164  H3K4me1_chr1-150428     chr1_1  NA      NA      FALSE
            chr1    235531  3596796 chr1_internal_44032     chr1_44032      chr1_1  chr1_30199      FALSE
            chr1    235531  3596796 chr1_internal_45833     chr1_45833      chr1_44032      chr1_41346      FALSE
            chr1    235531  16359528        chr1_internal_48404     chr1_48404      chr1_45833      chr1_46409      FALSE
            chr1    235531  16359528        chr1_internal_50073     chr1_50073      chr1_48404      chr1_48365      FALSE
            chr1    235531  32538399        chr1_internal_51024     chr1_51024      chr1_50073      chr1_49823      FALSE

        """

        Logger.get_instance().info( "read_crd_file: reading file %s" % ( self.crdFile) )

        if ".gz" in self.crdFile:  data = gzip.open(self.crdFile, "rt")
        else: data = open(self.crdFile, "r")
        
        crdInfo = {} # key -> crd_id, val -> list of children
        crdChro = {} # key -> crd_id, val -> chromosome
        count = 0
        _ = data.readline()
        for line in data.readlines():
  
            count+=1
            if count % 10000 == 0:
                Logger.get_instance().info( "read_crd_file: read %s entries" % ( count))

            line = line.strip()
            spl = line.split(" ")
 
            crdID = spl[3]
            if "internal" in crdID:
                crdID = crdID.split("_")
                chro = crdID[0]
                crdID = crdID[2]
                crdChro[crdID] = chro
                
            child1 = spl[1]
            child2 = spl[2]
            
            if child1 == "NA" and child2 == "NA":
                # peak entry
                child1 = int(spl[16])
                child2 = int(spl[17])
                crdID = spl[0]
            
            crdInfo[crdID] = [child1,child2]

        Logger.get_instance().info( "read_crd_file: read %s entries" % ( len(crdInfo) ) )
# 
        self.crdInfo = crdInfo
        self.crdChro = crdChro
        

    def read_gene_crd_link_file(self):
        """
        
        Read CRD-gene link file (format from Olivier Delaneau)
    
        
        Example format (space-separated, no header):
            ENSG00000000457.9 chr1_internal_36954 1.84551e-31 0.599998
            ENSG00000000457.9 chr1_internal_38026 3.24011e-05 0.234746
            ENSG00000000460.12 chr1_internal_38026 0.000373834 0.198244
            ENSG00000000938.8 chr1_internal_43387 4.27391e-25 0.540218
        
        """

        outFile = open(self.outFile, "w")

        Logger.get_instance().info( "read_gene_crd_link_file: reading file %s" % ( self.geneLinkFile) )

        if ".gz" in self.geneLinkFile:  data = gzip.open(self.geneLinkFile, "rt")
        else: data = open(self.geneLinkFile, "r")
         
        count = 0
        for line in data.readlines():
            count+=1
            if count % 1000 == 0:
                Logger.get_instance().info( "read_gene_crd_link_file: read %s entries" % ( count))
 
            line = line.strip()
            spl = line.split(" ")
            
            gene = spl[0]
            crdID = spl[7].split("_")[2]
                        
            # Search all child, store coordinates
            q = queue.Queue()
            q.put([crdID])
            allPeaks = []
            while not q.empty():
                n = q.get()
                if isinstance(n[0],int) and isinstance(n[1],int):
                    allPeaks.append(n)
                for s in n:
                    if s in self.crdInfo:
                        q.put(self.crdInfo[s])
        
            for peak in allPeaks:
                outFile.write("%s\t%s\t%s\t%s\n" % (self.crdChro[crdID], peak[0], peak[1], gene.split(".")[0]))

        outFile.close()

                
    def run(self):
        """
        Run functions in order
        """

        Timer.get_instance().step( "Read CRD file.." )        
        self.read_crd_file()

        Timer.get_instance().step( "Read link between genes and CRDs.." )        
        self.read_gene_crd_link_file()
        

if __name__ == "__main__":

    try:
    
        # Start chrono
        print ("STARTING " + SCRIPT_NAME)

        #===============================================================================
        # Get input arguments
        #===============================================================================
        parser = argparse.ArgumentParser(description= DESC_COMMENT) 
    
        # positional args
        parser.add_argument('crdFile', metavar='crdFile', type=str,
                             help='File with CRDs (tree).')
        parser.add_argument('geneLinkFile', metavar='geneLinkFile', type=str,
                             help='File with link between gene and CRDs.')
        parser.add_argument('outputFile', metavar='outputFile', type=str,
                             help='File where output will be written.')
                
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # Initialise class    
        run = ParseCRDBlueprint( args.crdFile, args.geneLinkFile, args.outputFile)

        run.run()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except Exception as e:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + str(e))

