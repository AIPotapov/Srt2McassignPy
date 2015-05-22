# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 22:31:26 2015

@author: Alexey
Using bmrb.py to read and process NMR-STAR files

Use bmrb module for reading NMR-STAR files
followed by processing into a Sequence object

14may made changes to accomodate a change in Sequence class API

My objects are in fact mutable,
 I haven't tried to make them immutable, but one just needs to be
 careful
 I have accidentelly discovered this when I called init_signal_table 
 appended the object, the called init... second time with different
 argument and appended the second  object to an array. It turned out that
 The first object in the table got changed
 I have to be careful with this and send copies of objects
 or alternatively create new object each time
"""


import sys
sys.path.append("../sans_python/")
import bmrb
import residue_test

# CONSTANTS
RESI_TYPES = ('ASP', 'THR', 'SER', 'GLU', 'PRO', 'GLY',
              'ALA', 'CYS', 'VAL', 'MET', 'ILE', 'LEU',
              'TYR', 'PHE', 'HIS', 'LYS', 'ARG', 'TRP',
               'GLN', 'ASN')
RESI_TYPES_BRIEF = ('D', 'T', 'S', 'E', 'P', 'G',
                    'A', 'C', 'V', 'M', 'I', 'L',
                    'Y', 'F', 'H', 'K', 'R', 'W',
                    'Q', 'N')
# end CONSTANTS
                    
def knock_out(seq, threshold=0.1, adj='no'):
    import random as rnd
    sample_size = int(threshold * len(seq))
    if adj is 'no':
        select_id = rnd.sample(seq.get_id_list(),
                               len(seq) - sample_size)        
    elif adj is 'yes':
        id_list = seq.get_id_list()
        id_list.sort()
        r = rnd.choice(id_list)
        if r + sample_size < len(seq):
            select_id = id_list[0: r]
            select_id.extend(id_list[r + sample_size: len(seq)])            
        else:
            select_id = id_list[r - len(seq) + sample_size : r]
    return seq.get_slice(select_id)

seq = residue_test.Sequence()
entry = bmrb.entry.fromFile('15156.srt')
#entry.printTree()

saveframeCS = entry.getSaveframesByCategory('assigned_chemical_shifts')
table_list = []
for saveframe in saveframeCS:
    loopCS = saveframe.getLoopByCategory('Atom_chem_shift')
    Data = loopCS.getDataByTag(['Seq_ID','Comp_ID','Atom_ID','Val'])
    for item in Data:
        seq.add_atom(int(item[0]), item[1], item[2], item[3])
    #print seq
    table = Table(seq)
    table.init_signal_table(table_type = 'N(CO)CACB')    
    table_list.append(table.copy())    
    table.init_signal_table(table_type = 'NCACX')    
    table_list.append(table.copy())
    #print table

#for item in table_list:
#    print item.print_MCASSIGN()
item = table_list[0]   
#table = table_list[0]
item.set_entry(0, resi_type = ['MET','ASP','GLU','GLY','C','A','X'])
item.set_entry(50, resi_type = ['MET','ASP','GLU','GLY','C','A','X'])
item.set_entry(100, resi_type = ['MET','ASP','GLU','GLY','C','A','X'])
#print item.print_mcassign()   
print seq
for i in range(100):
    seq2 = knock_out(seq, adj = 'yes', threshold = 0.2)
    print len(seq2, )
    #print seq2.get_id_list()
    #    print seq.get_id_list()

"""
# EARLY TESTS

ent15000 = bmrb.entry.fromDatabase(15000)
ent15000.printTree()

ent15000['entry_information']
ent15000['entry_information']['_Entry_author']

ent15000_2 = bmrb.entry.fromDatabase(15000)
del ent15000_2['entry_information']

bmrb.diff(ent15000, ent15000_2)

ent15000['entry_information']['_Entry_author'].columns
ent15000['entry_information']['_Entry_author'].data

a = ent15000['assigned_chem_shift_list_1']['_Atom_chem_shift'].data
a1 = ent15000['assigned_chem_shift_list_1']['_Chem_shift_experiment'].data
"""

