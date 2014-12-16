#!/usr/local/bin/python3

import sys
from LocusDB import Locus
from readsodict import readsodict 
import re
import sqlite3
# cmd line options: test500.vcf out.csv -a 1 -f 0.1

#####################################################
#####################################################
import argparse
parser = argparse.ArgumentParser('adds a vcf file to an SQLite database')
parser.add_argument("inFile", 
                    help="input file (.vcf) with gene annotation")
#parser.add_argument("outFile", 
#                    help="output file (.csv); for stdout type '-' ")

parser.add_argument("-d", "--dbPath", default = 'vcf.db',
                    help="a path to the database")
                    
parser.add_argument("-s", "--soterms", default="SO_terms.csv",
                    help="a path to a .csv file with SO terms and prior definition")

parser.add_argument("-f", "--frequencyFilter", type=float, default=0.0,
                    help="lower threshold of frequency to output")

parser.add_argument("-a", "--altCountFilter", type=int, default=0,
                    help="output reads with alternative counts strictly more than the given value")

parser.add_argument("-t", "--totCountFilter", type=int, default=0,
                    help="output reads with total counts strictly more than the given value")            
   
parser.add_argument("-u", "--sumAltCountFilter", type=int, default=0,
                    help="output loci with alternative counts strictly more than the given value")

parser.add_argument("-c", "--csvseparator", type=str, default= r';',
                    help="separator for the output .csv file [;]")     
                    
parser.add_argument("-l", "--label", type=str, default= r'vcf',
                    help="label of the database ['vcf']")     
                    
args = parser.parse_args()
#####################################################
def applyFilter(args, loc):
    flag = False
    totAlt = 0
    for ii in range(0, len(loc.pop)):
        totAlt += loc.pop[ii].altCount
        flag = flag or ((loc.pop[ii].altFrequency > args.frequencyFilter) \
            and (loc.pop[ii].altCount > args.altCountFilter) \
            and (loc.pop[ii].totCount > args.totCountFilter)) 
    flag = flag and totAlt > args.sumAltCountFilter
    # if flag:
    #     print('##########' , file=sys.stderr)
    return flag

###############################################################################

def initializeTables(conn, dbName, rowSqlTypes, countSqlTypes, chrNumber = 5):
     with conn:
        curs = conn.cursor()
        #######################################################################
        ## create the database tables (for each chromosome)
        for cc in range(1, chrNumber+1):
            rewriteTableChr(curs, dbName + '_vcf', rowSqlTypes, cc)
            rewriteTableChr(curs, dbName + '_mtCounts', countSqlTypes, cc)
            rewriteTableChr(curs, dbName + '_wtCounts', countSqlTypes, cc)

def rewriteTableChr(curs, tableName, rowSqlTypes, cc):
    query = 'DROP TABLE IF EXISTS %s_%u' % (tableName, cc)
    curs.execute(query)            
    
    query = ' CREATE TABLE %s_%u (' % (tableName, cc) + rowSqlTypes + ')'            
    curs.execute(query)           
###############################################################################
if len(args.soterms)>0 and not (args.soterms=='0') :
    SO_DICTIONARY = readsodict(args.soterms)
else:
    SO_DICTIONARY = []

#if not args.outFile == r'-':
#    sys.stdout = open(args.outFile, 'w')

f = open(args.inFile)

# columnnames = f.readline();

commre = re.compile(r"^[ ]*#.*");

loc = Locus("", SO_DICTIONARY);

###############################################################################
 
if (args.frequencyFilter and not args.altCountFilter):
    for line in f:
        if not (commre.match(line)):
            loc = Locus(line, SO_DICTIONARY);
            loc.printFields(args.csvseparator)
else:
    totLines = 0
    skippedLines = 0
    for line in f:
        if not (commre.match(line)):
            header = lastheaderline.split("\t");
            print("header: \n %s" % lastheaderline, file=sys.stderr, end='')
            
            numOfSamples = 0;
            repeatInfoFlag = False
            for sampleName in header[9:]:
                repeatInfoFlag = repeatInfoFlag or (sampleName.strip() == r'Repeat_Info')
                if not repeatInfoFlag:
                    print("sample no %u " % numOfSamples, ":'"+sampleName.strip()+"'", sep = None, file=sys.stderr)
                    numOfSamples += 1
            print('number of samples: %u' % numOfSamples , file=sys.stderr)
            break
        else:
            lastheaderline = line
        
    loc.numOfSamples = numOfSamples

    ###########################################################################
    ##  initialize the tables    

    print("writing to the database '%s'" % args.dbPath, file=sys.stderr)
    conn = sqlite3.connect(args.dbPath)
    dbName = args.label
    chrNumber = 5
    
    rowSqlTypes  = loc.sqlType()
    countSqlTypes = loc.pop[0].sqlType()
    initializeTables(conn, dbName, ', '.join(rowSqlTypes), ', '.join(countSqlTypes), 5)
    
    ###########################################################################    
    ## fill in the tables
    chrom = '0'
    with conn:
        curs = conn.cursor()
        print('-----------------------------', file=sys.stderr)
        print('-  filling in the database  -', file=sys.stderr)
        vcfQMarks = ','.join(['?'] * len(rowSqlTypes))
        countQMarks = ','.join(['?'] * len(countSqlTypes))        
        ## load vcf data into the database
        for line in f: 
                totLines += 1
                loc = Locus(line, SO_DICTIONARY, numOfSamples);
                if (applyFilter(args, loc)):
                    vcfRow = loc.getRowSqlite()
                    if (len(vcfRow[0])==1):   # ignore chloroplasts and mitochondria
                        query = u'INSERT INTO %s_vcf_%s VALUES ( %s )' % (dbName, vcfRow[0], vcfQMarks )
                        curs.execute(query, tuple(vcfRow[1]))                    
                        for ii in range(loc.numOfSamples):
                            countRow = loc.pop[ii].getRowSqlite()
                            query = u'INSERT INTO %s_%sCounts_%s VALUES ( %s )' % (dbName, loc.pop[ii].name, countRow[0], countQMarks)
                            curs.execute(query, tuple(countRow[1]) )
                        if not chrom == vcfRow[0]:
                            print('processing chromosome %s' % vcfRow[0], file=sys.stderr)
                        chrom = vcfRow[0]
                else:
                    skippedLines += 1
        print('-----------------------------', file=sys.stderr)        
        print( "out of %u lines," % totLines, "%u skipped" % skippedLines , "%u remained" % (totLines-skippedLines), file=sys.stderr)
###############################################################################


print('------------------', file=sys.stderr)  
with conn:
    curs = conn.cursor()
    curs.execute("SELECT * FROM sqlite_master WHERE type = 'table';")
    descr = curs.fetchall()
    r = curs.fetchone()
    for d in descr:         print(d[1])

f.close()

conn.close()

#if not args.outFile == r'-':
#    sys.stdout.close()
#    
sys.stdout = sys.__stdout__
