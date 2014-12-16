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
parser.add_argument("inFile", nargs='?', type=str, default='C24_snp_annotation_TAIR10.txt',
                   help="input file (.csv) with the gene IDs given in the 'geneID' column")
                   
parser.add_argument("outFile", nargs='?', type=str, default='vcf.db',
                    help="output file (.csv); for stdout type '-' ")

parser.add_argument("-l", "--label", type=str, default= r'C24',
                    help="label of the database ['C24']")     
                    
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
snpDB = r'call_method_75/call_method_75_TAIR9.csv'
dbName = args.label

###############################################################################
            
dbPath = 'vcf.db'
conn = sqlite3.connect(dbPath)
tt = oneecotypesnp(conn, args.label)

snpInputDB = args.inFile
tt.initializeSnpTable()

tt.populateSnpTable(snpInputDB, '\t')
    
print('------------------', file=sys.stderr)  
dprev = '#'
count = []
with conn:
    curs = conn.cursor()
    curs.execute("SELECT * FROM sqlite_master WHERE type = 'table';")
    descr = curs.fetchall()
    r = curs.fetchone()
    for d in descr:         
        count.append(d[1][-1:])
        d_n = d[1].split('_')[0] 
        if not (d_n== dprev):
            print(d_n, file=sys.stderr)
        else:
            # print(count, file=sys.stderr)
            count = []
        dprev = d_n


exmpl = tt.selectPosition(1, 5149, subtable = '', fields = 'ref, snp')

conn.close()
###############################################################################
