# -*- coding: utf-8 -*-
"""
Created on Mon May 25 16:28:29 2015

@author: Alexey

Generates the input files for MCASSIGN program
runs mcassign afterwards

Currently it works only for NCACX and N(CO)CACB files
but it can be extended later on
asfasfasfsf
"""
import os
import sys
sys.path.append("../sans_python/")
import bmrb
sys.path.append("./")
import residue_test

def write_connections():
    """
    generates connections
    currently the connections format is fixed for NCACB and N(CO)CACB tables
    but it may later be extended
    """
    s = ''
    s = '4\n'
    s = s + '1,2,4,4,0\n'
    s = s + '1,2,3,3,0\n'
    s = s + '1,2,2,2,0\n'
    s = s + '1,2,1,1,-1\n'
    return s
# end generate_connections

def write_AAsequence(sequence):
    """
    returns a string with amino acid characters
    """
    return sequence.print_seqAA()

def read_sequence_srt(filename):
    """
    parses the srt file using bmrb parser and creates a Sequence and Table
    data structures based on NCACX and N(CO)CACX types of experiments
    """
    seq = residue_test.Sequence()
    entry = bmrb.entry.fromFile(filename)
    # entry.printTree()

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
    
    # print seq.get_seqAA()
    # process the chemical shift entries
    saveframeCS = entry.getSaveframesByCategory('assigned_chemical_shifts')
    table_list = [] 
    # it turned out that there could be more than one 
    # assigned_chemical_shifts dataframe - 
    # it is not quite yet clear how I could choose the one i need
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
    return {'sequence': seq, 'tables' : 
             {'NCACX' : table_list[1], 'N(CO)CACB' : table_list[0]}}
# end read_sequence
             
def write_input(params):
    s = params['protein_sequence'] + '\n'
    s = s + str(params['number_of_signal_tables']) + '\n'
    s = s + str(params['signal_tables']['NCACX']) + '\n'
    s = s + str(params['signal_tables']['N(CO)CACB']) + '\n'
    s = s + str(params['connection_table']) + '\n'
    s = s + str(params['random_seed']) + '\n'
    s = s + str(params['good_connections']['initial']) + ' '
    s = s + str(params['good_connections']['final']) + '\n'
    s = s + str(params['bad_connections']['initial']) + ' '
    s = s + str(params['bad_connections']['final']) + '\n'
    s = s + str(params['edge_penalty']['initial']) + ' '
    s = s + str(params['edge_penalty']['final']) + '\n'
    s = s + str(params['signal_usage']['initial']) + ' '
    s = s + str(params['signal_usage']['final']) + '\n'
    s = s + str(params['step_number']['weighting']) + ' '
    s = s + str(params['step_number']['Monte Carlo']) + '\n'
    s = s + str(params['runs_number']) + '\n'
    s = s + str(params['output_prefix']) + '\n'
    s = s + str(params['output_file']) + '\n'    
    s = s + str(params['initial_assignments']) + '\n'    
    return s
    
# end write_input




# TEST CLIENT GENERATING CONNECTIONS FILE, SEQUENCE FILE
# MCASSIGN SIGNAL FILES

from params_mcassign import *

filename = '15156.srt'
a = read_sequence_srt(filename)
a['sequence'].change_name('b_' + filename[0:5])
sequence_name = a['sequence'].name

params['protein_sequence'] = sequence_name + '_seq.txt'
params['connection_table'] = sequence_name + '_conn.txt'
params['signal_tables'] = {'NCACX' : sequence_name + '_NCACX.txt',
                           'N(CO)CACB' : sequence_name + '_NCOCACB.txt'}
params['connection_table'] = sequence_name + '_conn.txt'
params['output_prefix'] = sequence_name[0:3] + 'dd'
params['output_file'] = sequence_name[0:3] + '.out'


print "*******SEQUENCE****************"
print write_AAsequence(a['sequence'])
with open(params['protein_sequence'], 'w') as f:
    f.write(write_AAsequence(a['sequence']))
    file_list = []
    file_list.append(params['protein_sequence'])

print "*******CONNECTIONS*************"
with open(params['connection_table'], 'w') as f:
    f.write(write_connections())
    file_list.append(params['connection_table'])

print "*******NCACX SIGNAL TABLE******"
with open(params['signal_tables']['NCACX'], 'w') as f:
    f.write(a['tables']['NCACX'].print_mcassign())
    file_list.append(params['signal_tables']['NCACX'])    

print "*******N(CO)CACB SIGNAL TABLE**"
with open(params['signal_tables']['N(CO)CACB'], 'w') as f:
    f.write(a['tables']['N(CO)CACB'].print_mcassign())
    file_list.append(params['signal_tables']['N(CO)CACB'])    


with open(sequence_name + '.in', 'w') as f:
    f.write(write_input(params))
    file_list.append(sequence_name + '.in')
    
file_list.append(params['output_file'])
    
    
#   BOOKKEEPING AFTER CALCULATION IS COMPLETE
import glob
file_list.extend(glob.glob(params['output_prefix'] + '*'))    

try:
    os.makedirs(a['sequence'].name)
except OSError:
    pass

import shutil
for file in file_list:
    try:
        shutil.move(file, a['sequence'].name + '/' + file)
    except IOError:
        print file + " does not exist"


