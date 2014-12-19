#!/usr/local/bin/python3

from vcfdb import GenomeVariantDbWriter
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
                    
parser.add_argument("-s", "--soPath", default="SO_terms.csv",
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
                    
parser.add_argument("-l", "--dbName", type=str, default= r'vcf',
                    help="label of the database ['vcf']")     
                    
args = parser.parse_args()
###############################################################################
    
# db = GenomeVariantDb( args.inFile, args.soPath, args.dbPath, args.dbName, chrNumber = 5)
db = GenomeVariantDbWriter( args = args, chrNumber = 5)

#if not args.outFile == r'-':
#    sys.stdout.close()
#    
# sys.stdout = sys.__stdout__
