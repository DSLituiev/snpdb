#!/usr/local/bin/python3.4
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:25:09 2014

@author: dima
"""
import sys
import sqlite3
from snpdbs import *

#####################################################
#####################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("inFile", 
                   help="input file (.csv) with the gene IDs given in the 'geneID' column")
                   
parser.add_argument("outFile", nargs='?', type=str, default='',
                    help="output file (.csv); for stdout type '-' ")

parser.add_argument("-l", "--label", type=str, default= r'vcf',
                    help="label of the database ['vcf']")     
                    
parser.add_argument("-r", "--rnaseq", default="./marc_rna_seq.txt",
                    help="a path to a table file with SO terms and prior definition")

parser.add_argument("-p", "--percentile", type=float, default=50.0,
                    help="threshold expressed as percentile of expression level")
                   
parser.add_argument("-c", "--columnsPositive", type=str, default='1,2',
                    help="column(s) of the rnaSeq data file, tissue/condition of interest")

parser.add_argument("-d", "--columnsNegative", type=str, default='3',
                    help="column(s) of the rnaSeq data file, negative control tissue/condition")
                    
parser.add_argument("-q", "--quantileflag", type=bool, default=False,
                    help="output reads with total counts strictly more than the given value")

parser.add_argument("-s", "--csvseparator", type=str, default= r';',
                    help="separator for the output .csv file [;]")
                    
parser.add_argument("-t", "--rnaSeparator", type=str, default= '\t',
                    help="separator for the input RNA expression table file ['\t']")
                    
parser.add_argument("-n", "--dataName", type=str, default= 'rnaFraction',
                    help="data name [rnaFraction]")
                    
parser.add_argument("-m", "--missingDataValue", type=float, default= float('NaN'),
                    help="value in case geneID is not present in the dataset [NaN]")                    
                    
args = parser.parse_args()

if not args.outFile == r'-':
    if not args.outFile == r'':
        sys.stdout = open(args.outFile, 'w')
    else:
        sys.stdout = open(args.inFile.replace('.csv','-ecotypeInfo.csv'), 'w')
        print('output file: %s' %args.inFile.replace('.csv','-ecotypeInfo.csv'),  file=sys.stderr) 


###############################################################################
currDir = r'/home/dima/scripts/snpdb/'
ecotypeInftsvoftsvilePath = r'call_method_75/call_method_75_info.tsv'
snpDB = r'call_method_75/call_method_75_TAIR9.csv'
dbName = args.label
###############################################################################
# conn = startSnpDatabase(currDir+ecotypeInftsvoftsvilePath, currDir+snpDB, dbPath = 'vcf.db')
conn = sqlite3.connect('vcf.db')
###############################################################################
# ecoT = tables(conn,'ecosnp')
# print('tables in the "vcf" DB: %s' % (dbName,  ecoT.showTables()),  file=sys.stderr) 
# print('number of SNPs in the ecoT DB: %u' % (dbName,  ecoT.countRows()),  file=sys.stderr) 

vcfT = tables(conn, dbName)

print('tables in the "vcf" DB: %s' % (dbName),  file=sys.stderr) 
print(vcfT.showTables(),  file=sys.stderr) 
vcfT.countRows('mtCounts')

c24T = oneecotypesnp(conn, 'C24')

with open(args.inFile) as f:
    header = next(f)
    chrInd = header.split(args.csvseparator).index('chr')
    posInd = header.split(args.csvseparator).index('pos')
    #refAlleleInd = header.split(args.csvseparator).index('refAllele')    
    # altAlleleInd = header.split(args.csvseparator).index('altAllele')
        
    endingStr = args.csvseparator+'\n'
#    print(header.rstrip('\n'),'', 'snpFrequencyInEcotypes', sep = args.csvseparator)    
    print(header.rstrip(';\n'),'', 'snpFrequencyInEcotypes', 'snpEcotypesInfo', sep = args.csvseparator)
    infoTable = ['','C24 or Col-0', 'most frequent allele (array)', 'polymorphic locus (array)', 'other than in Col']
    for line in f:   
        #start = time.clock()
        cols = line.split( ';') 
        #  print(line, file=sys.stderr)
        if not cols[posInd] == "":
            chromosome = int(cols[chrInd])
            pos = int(cols[posInd])                        
            (refAllele, altAllele) = vcfT.selectPosition( chromosome, pos,  'vcf', 'refAllele, altAllele')
            # c24Alleles = ecoT.selectPosition( chromosome, pos, 'refAllele, altAllele')
            
            c24Allele0 = c24T.selectPosition(chromosome, pos, subtable = '', fields = 'snp')

            c24Allele = getGenotypes(conn, chromosome, pos, 'C24')
            col0Allele = getGenotypes(conn, chromosome, pos, ecotype = 'Col-0')
            mfrAllele = getGenotypes(conn, chromosome, pos, ecotype = 'most frequent')
            mutFreq = getGenotypes(conn, chromosome, pos, 'all', altAllele)
#            print(chromosome, pos, '%3f' % mutFreq, refAllele, altAllele, col0Allele, c24Allele, mfrAllele, 
#            ''.join(getGenotypes(conn, chromosome, pos)),  file=sys.stderr)         
            if len(c24Allele0):
                # check c24 sequencing data
                if (altAllele == c24Allele0):
                    alleleFlag = 1
                else:
                    alleleFlag = 4
            elif len(c24Allele)>0:
                # check snp array data
                if (altAllele == c24Allele) or (altAllele == col0Allele)  :
                    alleleFlag = 1
                elif (refAllele == mfrAllele) and not (refAllele == col0Allele):
                    alleleFlag = 2
                else:
                    alleleFlag = 3
            else:
                alleleFlag = 0
                
#        if not alleleFlag:
#            print(line.rstrip('\n'),  mutFreq, sep = args.csvseparator, end = args.csvseparator + '\n')  
        if not alleleFlag == 1:          
            if not alleleFlag == 4:
                print(line.rstrip('\n'), mutFreq,  infoTable[alleleFlag], sep = args.csvseparator, end = args.csvseparator + '\n')
            else:
                print(line.rstrip('\n'), mutFreq,  'alt: %s|Col-0: %s| C24: %s' % (altAllele, refAllele, c24Allele0), sep = args.csvseparator, end = args.csvseparator + '\n')
                
print('finished successfully',  file=sys.stderr) 
        
if not args.outFile == r'-':
    sys.stdout.close()
    
sys.stdout = sys.__stdout__


###############################################################################
