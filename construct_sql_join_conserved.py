# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 15:27:46 2014

@author: dima
"""
def cons_table_chr(ss, chromosome):
    return 'Conservation_' + ss + '_%u'%chromosome


def cons_col(ss):
    return ss + '_cons'
    
species_list = ['Arabidopsis_lyrata' ,
'Aethionema_arabicum' ,
'Capsella_rubella' ,
'Brapa_197' ,
'Brassica_oleracea' ,
'Leavenworthia_alabamica' ,
'Sisymbrium_irio' ,
'Solanum_lycopersicum' ,
'Vitis_vinifera']

chromosome = 1;
# map(list,zip(species_list, numlist))

select_clause = '\n , '.join([ cons_table_chr(ss, chromosome) + '.conserved AS '+ cons_col(ss) for ss in species_list])

left_join_clause = '\n '.join([ ' LEFT JOIN ' + cons_table_chr(ss, chromosome) + ' ON '+ \
cons_table_chr(species_list[0], chromosome) + '.pos = ' + cons_table_chr(ss, chromosome) + '.pos ' \
for ss in species_list[1:] ])

where_clause = '\n OR '.join( [ '(' + cons_col(ss) + ' == 1 )' for ss in species_list] )

SQL_SENTENCE =  'SELECT ' + cons_table_chr(species_list[0], chromosome) + '.pos, ' + select_clause + \
'\n\n FROM ' + cons_table_chr(species_list[0], chromosome) + left_join_clause + \
'\n\n WHERE ' + where_clause

print(SQL_SENTENCE)