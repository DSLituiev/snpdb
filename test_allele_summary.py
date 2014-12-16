# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 15:34:51 2014

@author: dima
"""

from snpdbs import *

al_table = allele_summary()
al_table.calcFlags('A')
al_table.gts['most_freq_array'] = 'T'
al_table.gts['fukushima'] = 'A'
al_table.gts['c24_array'] = 'C'
al_table.calcFlags('A')
al_table.calcBitCode()
al_table.bitString

al_table.flags


new_table = allele_summary()
new_table.deciferBinary(al_table.bitString)

new_table.flags

new_table.printLegend()

shared_items = set(new_table.flags.items()) & set(al_table.flags.items())
assert ( len(shared_items)==len(al_table.flags) ) 