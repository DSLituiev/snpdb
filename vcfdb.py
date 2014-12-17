# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 14:29:27 2014

@author: dima
"""
from collections import defaultdict

def readsodict(soTermsFile):
    soDict = defaultdict(list)

    f = open(soTermsFile)
    
    for line in f:
        line = line.rstrip('\n')
        cols = line.split(';')
        cols[-1] = int(cols[-1])
        soDict[cols[0]].extend(cols[1:])
        
    f.close()
    return soDict
#####################################################################3
    
class GenomeVariantDb:
    
    def __init__(self, inFile, so_path, dbPath, dbName, chrNumber = 5):
        self.f = open(inFile)
        self.addSoDict(soPath)
        
        loc = Locus("", self.SO_DICTIONARY);
        self.numOfSamples = readVcfHeader( self.f )
        loc.numOfSamples = self.numOfSamples
        
        print("writing to the database '%s'" % dbPath, file=sys.stderr)
        self.conn = sqlite3.connect(dbPath)
        
        rowSqlTypes  = loc.sqlType()
        countSqlTypes = loc.pop[0].sqlType()
        
        initializeTables(self.conn, dbName, ', '.join(rowSqlTypes), ', '.join(countSqlTypes), 5)
        
    def addSoDict(self, soPath):
        if len(soPath)>0 and not (soPath=='0'):
            self.SO_DICTIONARY = readsodict(soPath)
        else:
            """ raise exception !"""
            self.SO_DICTIONARY = []
    
    def fill_tables(self):
        chrom = '0'
        with self.conn:
            curs = self.conn.cursor()
            print('-----------------------------', file=sys.stderr)
            print('-  filling in the database  -', file=sys.stderr)
            vcfQMarks = ','.join(['?'] * len(rowSqlTypes))
            countQMarks = ','.join(['?'] * len(countSqlTypes))        
            ## load vcf data into the database
            totLines = 0
            skippedLines = 0
            for line in f: 
                    totLines += 1
                    loc = Locus(line, SO_DICTIONARY, self.numOfSamples);
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
