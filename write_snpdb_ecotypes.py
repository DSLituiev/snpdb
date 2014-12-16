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
#parser.add_argument("inFile", 
#                   help="input file (.csv) with the gene IDs given in the 'geneID' column")
                   
parser.add_argument("outFile", nargs='?', type=str, default='vcf.db',
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
###############################################################################
currDir = r'/home/dima/scripts/snpdb/'
ecotypeInftsvoftsvilePath = r'call_method_75/call_method_75_info.tsv'
snpDB = r'call_method_75/call_method_75_TAIR9.csv'
dbPath = args.outFile
###############################################################################
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' , file=sys.stderr)

conn = startSnpDatabase(currDir+ecotypeInftsvoftsvilePath, currDir+snpDB, args.outFile)

with conn:
    query = "pragma integrity_check"
    curs = conn.cursor()
    curs.execute(query)
    a = curs.fetchall()
    print('test: %s' % query, file=sys.stderr)
    print(a[0][0], file=sys.stderr)
    

conn.commit()    
conn.close()



###############################################################################
