# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 14:18:38 2014

@author: Dmytro Lituiev
"""
import re
import sys

def strtofloat(s):
    if s == '':
        return 0
    else:
        return float(s)

class Pieces(object):
    """
    a basis class for `Locus` and `Segregants` classes. 
    Provides functionality for database handling.
    """
    def _getRowNamesSqlite_(self):
        ii = 0
        outList = [0]* ( len(self.outFields) )
        for ff in self.outFields:
            outList[ii] = ff
            ii += 1
        return outList[1:]
    
    def getRowSqlite(self):
        ii = 0
        outList = [0]* ( len(self.outFields) )
        for ff in self.outFields:
            outList[ii] = getattr(self, ff)
            ii += 1
        return (outList[0], outList[1:])
    
    def sqlType(self):
        fields = self._getRowNamesSqlite_()
        typeDict = {'int': 'INT', 'str': 'TEXT', 'float': 'REAL', 'none': 'NULL' }
        out = []
        for ff in fields:
            ftype = getattr(self, ff).__class__.__name__
            if ff == 'pos':
                colName = '%s %s PRIMARY KEY' % (ff, typeDict[ftype])
            else:
                colName = '%s %s' % (ff, typeDict[ftype])
            out.append(colName)
        return out
    

    
class Segregants(Pieces):    
    outFields = ['chr', 'pos', 'totCount', 'refCount', 'altCount']
    def __init__(self, name, cc, pos, refCount, altCount, refQuality = 0, altQuality = 0):        
        if not (refCount is None):            
            self.chr = cc
            self.pos = pos
            self.name = name
            self.refCount = refCount
            self.altCount = altCount
            self.totCount = refCount + altCount            
            self.refQuality =  refQuality    # Quality reference
            self.altQuality = altQuality # Quality alternative
            if self.totCount:
                self.altFrequency = float(self.altCount) / float(self.totCount)
            else:
                self.altFrequency = 0.0                    
        else:
            self.chr = 0
            self.pos = 0
            self.name = 'mt'
            self.refCount = 0
            self.altCount = 0
            self.totCount = 0
            self.altFrequency = 0

class Locus(Pieces):
    outFields = []
    segrPops = ["mt", "wt"]
    #===========
    # refCount = [0, 0];
    # altCount = [0, 0];
    # totCount = [0, 0];
    # refQuality = [0, 0];
    # altQuality = [0, 0];
    # altFrequency = [0.0, 0.0];
    #===========
    
    # mFlag = True; # True if the read comes from the mutant pool
    

        
    def __init__(self, line, SO_DICTIONARY, numOfSamples = 1, flag = 0):
        """
        INPUT:

        - `line` takes a line from a `.vcf` file
        - `SO_DICTIONARY` (read from a `.csv` file) provides annotation ranking,
        based on which an annotation to be stored in the database is selected
        - `numOfSamples` -- number of samples (one `vcf` file may carry information 
        about several samples). The sample-specific information is stored in the 
        `pop` attribute (of `Segregants` class)
        """
        self.numOfSamples = numOfSamples
        self.pop = [None]*numOfSamples
        #===========
         
        if len(SO_DICTIONARY)>0:
            self.outFields = ['chr', 'pos', 'rating', 'refAllele', 'altAllele',\
            'geneID', 'geneSO', \
            # , 'mutType',\
            'mutCDS', 'mutProt', 'mutPosCDS', 'mutPosProt', 'distance']    
        else:
            self.outFields = ['chr', 'pos', 'quality', 'refAllele', 'altAllele']
            
        if not line:  
            self.chr = 0;    
            self.pos = 0;
            self.quality = 0;
            self.refAllele = "";
            self.altAllele = "";        
            self.pop[0] = Segregants('mt', 0, 0, None, None)
            if self.numOfSamples>1:
                self.pop[1] = Segregants('wt', 0, 0, None, None)
            if len(SO_DICTIONARY)>0:
                self.geneID = "";
                self.geneSO = "";
                self.mutType = "";
                self.mutCDS = "";
                self.mutProt = "";
                self.mutPosCDS = 0;
                self.mutPosProt = 0;
                self.repeatName = "";
                self.repeatType = "";
                self.rating = 0;
                self.distance = 0.0# float("inf");
            
            return
        
        cols = line.split('\t')
        if (cols[0][0:3] == "Chr") or (cols[0][0:3] == "Chr"):
            self.chr = cols[0][3:]
        else:
            self.chr = cols[0]  # chromosome
        
        self.pos =   int(cols[1])    # Position
        self.refAllele = cols[3]     # ref allele    
        self.altAllele = cols[4]     # alt 
        
        if (cols[5] != "."):
            self.quality = cols[5]
        else:
            self.quality = 0
        
        ##
        def readSampleCounts(self, cols, ii):
            SUMMARY = cols[9+ii].split(":");
            # GT:DP:RO:QR:AO:QA:GL
            #  0  1  2  3  4  5  6
            # print('sample: %u'% ii, 'data: %s' % cols[9+ii], 'datalen %u' % len(cols[9+ii]), file=sys.stderr)
            if not (cols[9+ii]=='.') and len(cols[9+ii])>2:
                refCount = int(SUMMARY[2])     # Observations reference
                altCount = int(SUMMARY[4])     # Observations alternative
                # totCount = int(SUMMARY[1])     # Depth
                refQuality = SUMMARY[3]      # Quality reference
                altQuality = SUMMARY[5]      # Quality alternative
                return (refCount, altCount, refQuality, altQuality)
            else:
                return (0, 0, 0, 0)                
        
        
        if (flag==0 ):
            ii = 0
            while ii<self.numOfSamples:            
                (refCount, altCount, refQuality, altQuality) = readSampleCounts(self, cols, ii)
                self.pop[ii] = Segregants(self.segrPops[ii], self.chr, self.pos, refCount, altCount, refQuality, altQuality)
                ii +=1
        else:
            self.pop[0] = Segregants('mt', self.chr, self.pos, refCount, altCount, refQuality, altQuality)
            
        ## gene annotation column
        if len(SO_DICTIONARY)==0:
            return
        re1 = re.compile(";CSQ=.")
        try:
                INFO = re1.split(cols[7])[1];
        except:
                print(" probably there is no annotation in this '.vcf' file, or its format is wrong."\
                "\n please use Ensemble VEP to add gene annotation\n______________________________ ", file=sys.stderr)
                raise
        genes = INFO.split(',')
        ratingHit = 0
        geneHit = ['0']; # genes[0];
        for g in genes:
            gd = g.split("|", 12)
            if (len(gd)>1) :
                distance = float("inf");
                descr = gd[4].split('&');
                for dd in descr:
                    # print(dd, end = '\t')
                    rating = SO_DICTIONARY.get(dd);
                    # print(rating, end = '\t')
                    if (rating!= None) and ((rating[-1] > ratingHit) or ( (rating[-1] == ratingHit) and ( strtofloat(gd[-2]) < distance) )):
                        ratingHit = rating[-1];
                        gd[0] = SO_DICTIONARY.get(dd)[1]
                        geneHit = gd;
                        distance = strtofloat(gd[-2])
        if len(geneHit)<10:
            print( "geneHit entries are missing %s, %u" % (self.chr, self.pos) , file=sys.stderr)
            for i in range(len(geneHit)-1, 9): geneHit.append("")        


        self.geneSO     = geneHit[0][3:] # 'SO:' prefix is removed
        self.geneID     = geneHit[1]
        self.mutType    = geneHit[4]
        self.mutPosCDS  = geneHit[6]
        self.mutPosProt = geneHit[7]
        self.mutCDS     = geneHit[8]
        self.mutProt    = geneHit[9]
        self.distance   = geneHit[-2]
        self.rating     = ratingHit;
        
        ## repeat info
        if (len(cols)> 9+self.numOfSamples):
            REPEAT_INFO = cols[10+self.numOfSamples-1].split(";");
            self.repeatType = REPEAT_INFO[1];
            if  not (REPEAT_INFO[1] == "NO") :
                self.repeatName = REPEAT_INFO[0] # Repeat name
        
    def printFields(self, SEP):
        for name in self.outFields:
            val = getattr(self, name)
            print(val, end = SEP)
        for ind in range(0, self.numOfSamples):
            for name in self.outCountFields:
                val = getattr(self, name)
                print(val[ind], end = SEP)
                
        print('')
    
    def printFieldNames(self, SEP):
        print('num of samples: %u' % self.numOfSamples, file=sys.stderr)
        for name in self.outFields:
            print(name, end = SEP)
        for name in self.outCountFields:
            print('mt_'+name, end = SEP)
        if (self.numOfSamples > 1):
            for name in self.outCountFields:
               print('wt_'+name, end = SEP)
        
        print('')
        
#    def getRowSqlite(self):
#        ii = 0
#        outList = [0]* ( len(self.outFields) - 1)
#        for name in self.outFields[1:]:
#            outList[ii] = getattr(self, name)
#            ii += 1
#        return (self.chr, outList)
#            
#    def _getRowNamesSqlite_(self):
#        ii = 0
#        outList = [0]* ( len(self.outFields) -1)
#        
#        for name in self.outFields[1:]:
#            outList[ii] = name
#            ii += 1               
#        return outList
