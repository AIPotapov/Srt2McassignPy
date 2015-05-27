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
 
 22 may added knock_out function
 added reading in of the sequence from bmrb file (it should be later
 used to construct sequence file of for the mcassign)
 the dictionary with sequence will be stored along with the Sequence
 
 25 may 
 1. added Mcassign class
 which data structure is constructed by parsing mcassign output file
 the data structure can be compared to another data structure
 the results can be used to give a report about similarity of difference
 of assignments with varying length
 2. Move knock_out to Sequence method
 
"""


import sys
sys.path.append("../sans_python/")
import bmrb
sys.path.append("./")
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
                    


# MAIN BODY OF THE PROGRAM
seq = residue_test.Sequence()
entry = bmrb.entry.fromFile('15156.srt')
entry.printTree()


# process the residue type entries to get full sequence
saveframeAA = entry.getSaveframesByCategory('entity')
for saveframe in saveframeAA:
    loopAA = saveframe.getLoopByCategory('Entity_comp_index')
    Data = loopAA.getDataByTag(['ID', 'Comp_ID'])    
    comp_ID = [row[1] for row in Data]
    ID = [row[0] for row in Data]
    numID = []
    for item in ID:
        numID.append(int(item))
    sequenceAA = dict(zip(numID, comp_ID))
  

seq.set_seqAA(sequenceAA)
print seq.get_seqAA()


# process the chemical shift entries
saveframeCS = entry.getSaveframesByCategory('assigned_chemical_shifts')
table_list = []

for saveframe in saveframeCS:
    loopCS = saveframe.getLoopByCategory('Atom_chem_shift')
    Data = loopCS.getDataByTag(['Seq_ID','Comp_ID','Atom_ID','Val'])
    for item in Data:
        seq.add_atom(int(item[0]), item[1], item[2], item[3])
    #print seq
    table = residue_test.Table(seq)
    table.init_signal_table(table_type = 'N(CO)CACB')    
    table_list.append(table.copy())    
    table.init_signal_table(table_type = 'NCACX')    
    table_list.append(table.copy())
    #print table





""" 
 this is just for testing
"""
#for item in table_list:
#    print item.print_MCASSIGN()
item = table_list[0]   
#table = table_list[0]
item.set_entry(0, resi_type = ['MET','ASP','GLU','GLY','C','A','X'])
item.set_entry(50, resi_type = ['MET','ASP','GLU','GLY','C','A','X'])
item.set_entry(100, resi_type = ['MET','ASP','GLU','GLY','C','A','X'])
#print item.print_mcassign()   
#print seq
for i in range(10):
    seq2 = seq.knock_out(adj = 'yes', threshold = 0.2)
    print len(seq2)
print item.print_seqAA()
print item.sequence.get_seqAA()

print "************TESTING MCASSIGN PARSING**********"

m1 = residue_test.Mcassign('HBVdd001.txt')
m2 = residue_test.Mcassign('HBVdd002.txt')
comp = m1.compare(m2)
print len(comp['similar'])
print len(comp['different'])


    #print seq2.get_id_list()
    #    print seq.get_id_list()