# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 13:11:16 2014

@author: dima
"""
import sys
import sqlite3
from bitmagic import *
# snpdb

###########################################
## initialize the snps table
def initializeSnpTable(conn, chrNumber = 5):
    
    dbases = [''] * chrNumber
    
    with conn:
        curs = conn.cursor()
        
        for cc in range(1,chrNumber+1):
            dbases[cc-1] = 'ecosnp_%u'% cc
            query = 'DROP TABLE IF EXISTS '+ dbases[cc-1]
            print(query, file=sys.stderr)
            curs.execute(query)
        
            query = ' CREATE TABLE ' + dbases[cc-1] + '''(
            pos   INT    PRIMARY KEY,
            genotypes TEXT
            )
            '''
            curs.execute(query)

###############################################################################
#def startSnpDatabase(conn, ecotypeInftsvoftsvilePath, snpDB, chrNumber = 5):
def startSnpDatabase(ecotypeInftsvoftsvilePath, snpDB, dbPath = 'snp_ecotype.db', chrNumber = 5):
    conn = sqlite3.connect(dbPath)

#    ###########################################     
    if not checkTable(conn, 'ecotypes'):        
        print('creating "ecotypes" table ', file=sys.stderr)
        # create an ecotype table
        (nn, ec) = readEcotypeDescription(ecotypeInftsvoftsvilePath)
        populateEcotypeTable(conn, ec, nn)    
    ###########################################     
    if not checkTable(conn, 'ecosnp_%u' % chrNumber ):
        ## initialize the snps table
        initializeSnpTable(conn, chrNumber)
#        ## populate the table
        populateSnpTable(conn, snpDB) 
        print('snp table is done', file=sys.stderr)
#    ###########################################     
#        a = getGenotypes(conn, 1, 1024, 'all', '')
#        print(a, file=sys.stderr)
        a = checkTable(conn, tablename = 'ecosnp_5')
        print('test: %u'% a, file=sys.stderr)
    return conn

###############################################################################
def readEcotypeDescription(ecotypeInftsvoftsvilePath):

    ftsv = open(ecotypeInftsvoftsvilePath)
    header = next(ftsv)
    # array_id_i = header.split('\t').index('array_id')
    ecotype_id_i = header.split('\t').index('ecotype_id')
    nativename_i = header.split('\t').index('nativename')
    
    nn = []
    ec = []
    
    for line in ftsv:   
        cols = line.rstrip('\n').split( '\t') 
        # ar.append( int(cols[array_id_i]) )
        ec.append( int(cols[ecotype_id_i]) )
        nn.append( cols[nativename_i] )
    
    # ecotypeInftsvo = dict(zip(ec, nn))
    # ecotypeIds =  dict(zip(nn, ec))
    
    print('ecotype description has been read successfully, %u entries, %u unique' \
    % ( len(ec), len(set(ec))), file=sys.stderr)
    
    ftsv.close()
    
    return (nn, ec)
## ecotype table
###########################################  
def checkTable(conn, tablename):    
    with conn:
        query = "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='%s'" % tablename
        curs = conn.cursor()
        curs.execute(query)
        return curs.fetchall()[0][0]
    
###########################################  
# create the ecotype table
def populateEcotypeTable(conn, ec, nn):
    with conn:           
        query = 'DROP TABLE IF EXISTS ecotypes'        
        curs = conn.cursor()
        curs.execute(query)
        
        query = '''
        CREATE TABLE ecotypes(
        id   INT    PRIMARY KEY,      
        arrayID   INT,
        et VARCHAR(16)    
        )
        '''
        curs.execute(query)
        
        query = 'INSERT INTO ecotypes VALUES (?,?,?)'
        for ii in range(1, len(nn)):
            curs.execute(query,  (ii, ec[ii], nn[ii]) )
###########################################
 
def populateSnpTable(conn, snpDB, separator = ','):
    
    fcsv = open(snpDB)
    next(fcsv)
    header = next(fcsv)
    hlen = len(header.split(separator))
    
    with conn:
        curs = conn.cursor()
        for line in fcsv:
            cols = line.rstrip('\n').split(separator);
            if hlen == len(cols):
                query = 'INSERT INTO ecosnp_%s VALUES (?,?)' % cols[0]
                b = ''.join(cols[2:])
                curs.execute(query,(cols[1], b))
            else:
                 print('the length of the row (%u) does not match the length of the header (%u)' \
                 % (len(cols), hlen), file=sys.stderr)
    
    fcsv.close()
###########################################

def getEcotypeArrayID(conn, eco0):
    with conn:
        curs = conn.cursor()
        curs.execute("SELECT arrayID FROM ecotypes WHERE et = '" + eco0+ "'")          
        return curs.fetchall()[0]

def getEcotypeID(conn, eco0):
    with conn:
        curs = conn.cursor()
        curs.execute("SELECT id FROM ecotypes WHERE et = '" + eco0+ "'")          
        return curs.fetchall()[0]        

###########################################        
from collections import Counter

def getGenotypes(conn, chrs, position, ecotype = 'all', alleleQuery = ''):
    # ecotype = u'C24'
    dbName = 'ecosnp'
    with conn:
        curs = conn.cursor()
        query = 'SELECT * FROM %s_%u WHERE pos = %u' %  (dbName, chrs, position)
        curs.execute(query)
        a = curs.fetchall()
        if len(a) == 0:
            if alleleQuery == '':
                return ''
            else:
                return 0
        genotypes = list(a[0][1])
        if alleleQuery == '':
            if ecotype == 'all':
                return genotypes
            elif ecotype == 'most frequent':
                stats = Counter(genotypes)
                return max(stats, key=stats.get)
            else:
                query = 'SELECT id FROM ecotypes WHERE et = ?'
                curs.execute(query, (ecotype,) )
                a = curs.fetchall()[0][0]
                return genotypes[a]
        else:
             stats = Counter(genotypes)
             return float(stats[alleleQuery]) / float(len(genotypes))
            
###############################################################################
class tables(object):
    def __init__(self, conn, nameBase):
        self.connection = conn
        self.nameBase = nameBase
#        try:
#            cursor = conn.cursor()        
#            cursor.execute("SELECT VERSION()")
#            results = cursor.fetchone()
#            # Check if anything at all is returned
#        except sqlite3.Error:
#            print("ERROR IN CONNECTION",  file=sys.stderr) 
        # self.ex = self.connection.cursor().execute

    def getname(self, subtable = ''):
        if subtable == '':
            self.name = self.nameBase
        else:
            self.name = '%s_%s' % (self.nameBase, subtable)

        return self.name
        
    def countRows(self, subtable = '', chrN = 5):
        a=0;
        print('number of SNPs in the DB: ',  file=sys.stderr) 
        for chromosome in range(1, chrN+1):
            query = 'SELECT COUNT(*) FROM %s' % (self.getnamechr( chromosome, subtable) )
            curs = self.connection.cursor()
            curs.execute(query)
            # a += curs.fetchall()
            print('chr %u: %u' % (chromosome, curs.fetchall()[0][0]),   file=sys.stderr) 
        # return a
        
        
    def showTables(self):
        with self.connection:
            query = "SELECT name FROM sqlite_master WHERE type='table';"
            curs = self.connection.cursor()
            curs.execute(query)
            a = curs.fetchall()
        return a
        
        
    def getnamechr(self, chromosome, subtable = ''):
        if isinstance(chromosome, int):
            chromosome = '%u' % chromosome
        
        self.getname(subtable)            
        self.nameChr = self.name + '_%s'%chromosome
        return self.nameChr
    
    def selectPosition(self, chromosome, pos, subtable = '', fields = '*'):
        if isinstance(pos, int):
            pos = '%u' % pos
        with self.connection:
            query = 'SELECT %s FROM %s WHERE pos = %s' % (fields, self.getnamechr( chromosome, subtable) , pos)
            curs = self.connection.cursor()
            try:
                curs.execute(query)                
                a = curs.fetchall()
            except:
                a = []
            
        if ',' in fields and ( (a == []) ):
            print('no one entry has been found, chr %u, position %s, returned: %s' % (chromosome, pos, a),  file=sys.stderr) 
            a = (0,) * len(fields.split(','))
            # print('returnung: ', end='', file=sys.stderr)
            # print( a, file=sys.stderr)
            return a
            
            # print('finished successfully',  file=sys.stderr) 
        
        if (a == []):
            return a        
        if not ((fields == '*') or ',' in fields):
             return a[0][0]
        else:
             return a[0]
###############################################################################
class oneecotypesnp(tables):
    ## class for C24 ecotype data
    def initializeSnpTable(self, chrNumber = 5):
        dbases = [''] * chrNumber
        
        with self.connection:
            curs = self.connection.cursor()
            
            for cc in range(1,chrNumber+1):
                dbases[cc-1] = '%s_%u'% (self.nameBase, cc)
                query = 'DROP TABLE IF EXISTS '+ dbases[cc-1]
                print(query,  file=sys.stderr) 
                curs.execute(query)
            
                query = ' CREATE TABLE ' + dbases[cc-1] + '''(
                pos   INT    PRIMARY KEY,
                ref CHAR,
                snp CHAR
                )
                '''
                curs.execute(query)
    ###########################################################################
    def populateSnpTable(self, snpInputDB, separator = '\t'):
        
        fcsv = open(snpInputDB)
        ## NO HEADER! 
#        next(fcsv)
#        header = next(fcsv)
#        hlen = len(header.split(separator))
        
        with self.connection:
            curs = conn.cursor()
            for line in fcsv:
                cols = line.rstrip('\n').split(separator);
#                if hlen == len(cols):
                query = 'INSERT INTO %s_%s VALUES (?,?,?)' % (self.nameBase, cols[1])
                curs.execute(query,(cols[2], cols[3], cols[4]))
#                else:
#                     print('the length of the row (%u) does not match the length of the header (%u)' \
#                     % (len(cols), hlen), file=sys.stderr)
        
        fcsv.close()

###############################################################################

class allele_summary:
    """a summary table of alleles; represents alleles as a binary code >> integer
    """
    ecotypes = ['fukushima', 'c24_genomic', 'c24_array', 'most_freq_array', 'polymorphic_array']

    def __init__(self):
        self.flags = {}
        self.gts = {}
        for ee in self.ecotypes:
            self.gts[ee] = ''
            self.bitString = int(0);
            
    def calcFlags(self, altAllele):
        for ee in self.ecotypes[0:-1]:
            if self.gts[ee] == altAllele:
                self.flags[ee] = True
            else:
                self.flags[ee] = False
            
        if len(self.gts['most_freq_array'])>0:
            self.flags['polymorphic_array'] = True
        else:
            self.flags['polymorphic_array'] = False
    
    def calcBitCode(self):
        self.bitString = int(0)
        for ll in range(0,len(self.ecotypes)):
            self.bitString = self.bitString + self.flags[self.ecotypes[ll]]*2**ll
        return self.bitString 

    def deciferBinary(self, num):
        assert (len(self.flags)==0), 'non-empty genotype flags!\n do you really want to over-write them?\n then creat a new instance'
        bv = get_bin_vect(int(num), len(self.ecotypes))

        for bb, ee in zip(bv, self.ecotypes):
            self.flags[ee] = bb
            
    def printLegend(self, file_pt=sys.stderr):
        print('Legend for the "snpEcotypesInfo" field:',  file=file_pt)
        num_eco = len(self.ecotypes)
        form_str = '\t'.join( ['%s']*(num_eco+1) )
        form_num = '\t'.join( ['%u']*(num_eco+1) )
        print( form_str % (('code',)+tuple(self.ecotypes)) ,  file=file_pt)
        
        for bb in range(0,2**num_eco):
            bv = get_bin_vect(int(bb), num_eco)
            print(form_num % ((bb,) + tuple(bv)),  file=file_pt) 
#########################################################