#!/usr/local/bin/python3.4
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:25:09 2014

@author: dima
"""
import sys
import sqlite3
from snpdbs import *

import os

try:
    CURR_DIR = os.path.dirname(os.path.realpath(__file__))
except NameError:
    CURR_DIR = '.'
DB_PATH = os.path.join(CURR_DIR, 'vcf.db')

#####################################################

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("inFile", 
                   help="input file (.csv) with the gene IDs given in the 'geneID' column")
                   
parser.add_argument("outFile", nargs='?', type=str, default='',
                    help="output file (.csv); for stdout type '-' ")

parser.add_argument("-s", "--csvseparator", type=str, default= r';',
                    help="separator for the output .csv file [;]")
                    
parser.add_argument("-m", "--missingDataValue", type=float, default= float('NaN'),
                    help="value in case geneID is not present in the dataset [NaN]")

parser.add_argument("-b", "--database", type=str, default = DB_PATH,
                    help="value in case geneID is not present in the dataset [NaN]")    
                    
args = parser.parse_args()

###############################################################################
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
###############################################################################
# conn = startSnpDatabase(currDir+ecotypeInftsvoftsvilePath, currDir+snpDB, dbPath = 'vcf.db')
if os.path.isfile(args.database):
    conn = sqlite3.connect(args.database)
else:
    raise OSError('database file not found')
###############################################################################
def forceInt(x):
    try:
        y = int(x)
    except:
        y = 0
    return y
###############################################################################
# ecoT = tables(conn,'ecosnp')
# print('tables in the "vcf" DB: %s' % (dbName,  ecoT.showTables()),  file=sys.stderr) 
# print('number of SNPs in the ecoT DB: %u' % (dbName,  ecoT.countRows()),  file=sys.stderr) 

c24T = oneecotypesnp(conn, 'C24')

fuku = oneecotypesnp(conn, 'MEA250n1_vcf')


with open(args.inFile) as f:
    header = next(f)    
    splitHeader = header.split(args.csvseparator)
    colInds = {}    
    for ii in range(0,len(splitHeader)):
        colInds[splitHeader[ii]] = ii        

    #refAlleleInd = header.split(args.csvseparator).index('refAllele')    
    # altAlleleInd = header.split(args.csvseparator).index('altAllele')
        
    endingStr = args.csvseparator+'\n'
#    print(header.rstrip('\n'),'', 'snpFrequencyInEcotypes', sep = args.csvseparator)    
    print(header.rstrip(';\n'), args.csvseparator, sep ='', end = '')
    print('snpFrequencyInEcotypes', 'snpEcotypesInfo', sep = args.csvseparator)
    # infoTable = ['','C24 or Col-0', 'most frequent allele (array)', 'polymorphic locus (array)', 'other than in Col']
    for line in f:
        mutFreq = 0 ;
        #start = time.clock()
        cols = line.split( ';') 
        #  print(line, file=sys.stderr)
        chromosome = forceInt(cols[colInds['chr']])
         
        ###### !!!!
        
        if not cols[colInds['pos']] == "" and not (chromosome== 0):            
            pos = int(cols[colInds['pos']])
            refAllele = str(cols[colInds['refAllele']])
            altAllele = str(cols[colInds['altAllele']])
            
            al_table = allele_summary()
            
            al_table.gts['c24_genomic'] = c24T.selectPosition(chromosome, pos, subtable = '', fields = 'snp')
            al_table.gts['fukushima']   = fuku.selectPosition(chromosome, pos, subtable = '', fields = 'altAllele')
            al_table.gts['c24_array']       = getGenotypes(conn, chromosome, pos, 'C24')
            al_table.gts['most_freq_array'] = getGenotypes(conn, chromosome, pos, ecotype = 'most frequent')
            mutFreq = getGenotypes(conn, chromosome, pos, 'all', altAllele)
#            print(chromosome, pos, '%3f' % mutFreq, refAllele, altAllele, col0Allele, c24Allele, mfrAllele, 
#            ''.join(getGenotypes(conn, chromosome, pos)),  file=sys.stderr)         

            al_table.calcFlags(altAllele)
            alleleFlag = al_table.calcBitCode()
 
                
#        if not alleleFlag:
#            print(line.rstrip('\n'),  mutFreq, sep = args.csvseparator, end = args.csvseparator + '\n')  
            print(line.rstrip(';\n'), args.csvseparator, sep ='', end = '')
            print(mutFreq,  alleleFlag, sep = args.csvseparator, end = args.csvseparator + '\n')
            
print('finished successfully',  file=sys.stderr) 

al_table.printLegend()
        
if not args.outFile == r'-':
    sys.stdout.close()
    
sys.stdout = sys.__stdout__


###############################################################################
