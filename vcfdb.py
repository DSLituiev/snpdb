# -*- coding: utf-8 -*-
"""
Data from a `vcf` file is put into a `sqlite` database with a following structure:

- Separate chromosomes go into separate tables;

- Each chromosome data is split into two types of tables:
  
  + Locus description table: 
    alt & ref letters, annotation
  
  + Sample description table[s]: 
    alt & ref counts for each sample
    
- All tables have nucleotide position as the key column

Created on Wed Dec 17 14:29:27 2014
@author: D S Lituiev
"""
from LocusDB import Locus
from readsodict import readsodict 
import sqlite3
import re
import sys
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
def readVcfHeader(f):
    commre = re.compile(r"^[ ]*#.*")
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
        
    return numOfSamples
###########################################################################
def applyFilter(Filter, loc):
    flag = False
    totAlt = 0
    for ii in range(0, len(loc.pop)):
        totAlt += loc.pop[ii].altCount
        flag = flag or ((loc.pop[ii].altFrequency > Filter['frequency'] ) \
            and (loc.pop[ii].altCount > Filter['altCount']) \
            and (loc.pop[ii].totCount > Filter['totCount'])) 
    flag = flag and totAlt > Filter['sumAltCount']
    # if flag:
    #     print('##########' , file=sys.stderr)
    return flag
###############################################################################
class GenomeVariantDbWriter:
    
    def __init__(self, inFile="", soPath="", dbPath="", dbName="", chrNumber = 5, Filter = {}, args=None):
        if not (args is None):
            self.inFile = args.inFile
            self.add_so_dict(args.soPath)
            self.dbPath = args.dbPath
            self.dbName = args.dbName           
            self.Filter = {}        
            self.Filter['frequency'] = args.frequencyFilter
            self.Filter['altCount'] = args.altCountFilter
            self.Filter['totCount'] = args.totCountFilter
            self.Filter['sumAltCount']= args.sumAltCountFilter
        else:
            self.inFile = args.inFile
            self.add_so_dict(soPath)
            self.Filter = Filter
            self.dbPath = dbPath
        
        self.init_db()
        self.fill_tables()
        self.check_tables()
        
    def init_db(self):
        loc = Locus("", self.SO_DICTIONARY);
        """ VCF file is open here ! """
        self.vcf = open(self.inFile, 'r')
        self.numOfSamples = readVcfHeader( self.vcf )
        loc.numOfSamples = self.numOfSamples
        # self._vcf_header_end_ = self.vcf.tell()
        # self.vcf.close()
        
        print("writing to the database '%s'" % self.dbPath, file=sys.stderr)
        
        self.rowSqlTypes  = loc.sqlType()
        self.countSqlTypes = loc.pop[0].sqlType()
        
        with sqlite3.connect(self.dbPath) as self.conn:
            initializeTables(self.conn, self.dbName, ', '.join(self.rowSqlTypes), ', '.join(self.countSqlTypes), 5)

        
    def add_so_dict(self, soPath):
        if len(soPath)>0 and not (soPath=='0'):
            self.SO_DICTIONARY = readsodict(soPath)
        else:
            """ raise exception !"""
            self.SO_DICTIONARY = []
    
    def fill_tables(self):
        chrom = '0'
        # self.vcf.open('r')
        # self.vcf.seek(self._vcf_header_end_)
        
        with sqlite3.connect(self.dbPath) as self.conn:
            curs = self.conn.cursor()
            print('-----------------------------', file=sys.stderr)
            print('-  filling in the database  -', file=sys.stderr)
            vcfQMarks = ','.join(['?'] * len(self.rowSqlTypes))
            countQMarks = ','.join(['?'] * len(self.countSqlTypes))        
            " load vcf data into the database "
            totLines = 0
            skippedLines = 0
            for line in self.vcf: 
                    totLines += 1
                    loc = Locus(line, self.SO_DICTIONARY, self.numOfSamples);
                    if (applyFilter(self.Filter, loc)):
                        chrom_vcf, vcfRow = loc.getRowSqlite()
                        if (len(chrom_vcf)==1):   # ignore chloroplasts and mitochondria
                            query = u'INSERT INTO %s_vcf_%s VALUES ( %s )' % (self.dbName, chrom_vcf, vcfQMarks )
                            curs.execute(query, tuple(vcfRow))                    
                            for ii in range(loc.numOfSamples):
                                countRow = loc.pop[ii].getRowSqlite()
                                query = u'INSERT INTO %s_%sCounts_%s VALUES ( %s )' % \
                                (self.dbName, loc.pop[ii].name, countRow[0], countQMarks)
                                curs.execute(query, tuple(countRow[1]) )
                            if not chrom == chrom_vcf:
                                print('processing chromosome %s' % chrom_vcf, file=sys.stderr)
                            chrom = chrom_vcf
                    else:
                        skippedLines += 1
            self.conn.commit()
            # curs.execute("COMMIT TRANSACTION");
            print('-----------------------------', file=sys.stderr)        
            print( "out of %u lines," % totLines, "%u skipped" % skippedLines , "%u remained" % (totLines-skippedLines), file=sys.stderr)
            
        self.vcf.close()

        
    def check_tables(self):
        print('------------------', file=sys.stderr)  
        with sqlite3.connect(self.dbPath) as self.conn:
            curs =  self.conn.cursor()
            curs.execute("SELECT * FROM sqlite_master WHERE type = 'table';")
            descr = curs.fetchall()
            r = curs.fetchone()
            for d in descr:         print(d[1])
            
###############################################################################

def readPosition(conn, dbName, chromosome, pos, field = '*'):
    with conn:
        query = 'SELECT %s FROM %s_%u WHERE pos = %u' % (field, dbName, chromosome, pos)
        curs = conn.cursor()
        curs.execute(query)        
        a = curs.fetchall()
        if not len( a ):
            return None
            
        if not field == '*':
            return a[0][0]
        else:
            return a[0]

class GenomeVariantDbReader:
    
    def __init__(self, soPath="", dbPath="", dbName="", chrNumber = 5, Filter = {}, args=None):
        if not (args is None):
            self.inFile = args.inFile
            self.add_so_dict(args.soPath)
            self.dbPath = args.dbPath
            self.dbName = args.dbName           
            self.Filter = {}        
        else:
            self.inFile = args.inFile
            self.add_so_dict(soPath)
            self.Filter = Filter
            self.dbPath = dbPath
            
    def read_counts(self, sample = "mt"):
        table_name = self.dbName + '_'+ sample + 'Counts'
        with sqlite3.connect(self.dbPath) as self.conn:
            readPosition(self.conn, table_name, chromosome, pos, field = '*')
        
