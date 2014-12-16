#!/usr/local/lib/python3.4
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 10:31:26 2014

@author: dima
"""

import sys
from LocusDB import Locus

sys.path.append('../vcftocsv_python/')

from readsodict import readsodict 
import re
import sqlite3
from snpdbs import *

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


###############################################################################
    

dbPath = 'vcf.db'
conn = sqlite3.connect(dbPath)
dbName = 'MO8'
print('------------------', file=sys.stderr)

table = tables(conn, dbName)

dbTables = []
with conn:
    curs = conn.cursor()
    curs.execute("SELECT * FROM sqlite_master WHERE type = 'table';")
    descr = curs.fetchall()
    r = curs.fetchone()
    for d in descr:  
        dbTables.append(d[1])
        print(d[1])
    
    chrTables = [ [] for i in range(5)]
    for ii in range(1,6):
        for t in dbTables:
            try:
                if (int(t[-1] ) == ii) and (t[0:len(dbName)] == dbName):
                    chrTables[ii-1].extend([t])
            except:
                pass
                    
    print(chrTables)
    
    curs.execute("select * from %s" % table.getnamechr(1, 'mtCounts') )
    colNames = [f[0] for f in curs.description]
    
    query = u'SELECT * FROM %s WHERE pos = %u' % (chrTables[0][0] , 1024)
    curs.execute(query)
    descr = curs.fetchall()
    

print('------------------', file=sys.stderr)