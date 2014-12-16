# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 14:21:00 2014

@author: Dmytro Lituiev
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