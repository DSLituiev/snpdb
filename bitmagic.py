# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 15:36:58 2014

@author: dima
"""

#####################################################
def get_bit(byteval,idx):
    return ((byteval&(1<<idx))!=0);

def get_bin_vect(byteval, bit_depth = 8):
    bit_list = [False]*bit_depth
    
    for n in range(0,bit_depth):
        bit_list[n] = get_bit( byteval, n)
    return bit_list
