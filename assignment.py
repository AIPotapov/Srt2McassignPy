# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 22:26:14 2015

Author: Alexey Potapov
Module contains classes for constructing MCASSIGN input files
and parsing MCASSIGN output files

some methods are not safe because no format or null
checking takes place

Progress: 
3 may added get_residue
and get_id_list methods
implementation of Table Class - not tested yet

4 may modify Residue class to include chem_shift_err
first working version of the Table class
it not yet thoroughly tested and outline may change in the future

14 may
update the implementation of _Entry class

25 may
move Mcassign and its helper class _TableLine from srt2mcassign

28 may 
moved all the procedures in run_mcassign.py program into a Run class
such arrangement could allow some free play with input sequences etc.
"""

import sys
import re
import bmrb
import os

# CONSTANTS
_RESI_TYPES = ('ASP', 'THR', 'SER', 'GLU', 'PRO', 'GLY',
               'ALA', 'CYS', 'VAL', 'MET', 'ILE', 'LEU',
               'TYR', 'PHE', 'HIS', 'LYS', 'ARG', 'TRP',
               'GLN', 'ASN')
_RESI_TYPES_BRIEF = ('D', 'T', 'S', 'E', 'P', 'G',
                     'A', 'C', 'V', 'M', 'I', 'L',
                     'Y', 'F', 'H', 'K', 'R', 'W',
                     'Q', 'N')
_RESI_DICT = dict(zip(_RESI_TYPES, _RESI_TYPES_BRIEF))
_RESI_DICT_REVERSE = dict(zip(_RESI_TYPES_BRIEF, _RESI_TYPES))
# end CONSTANTS


class Residue(object):
    """
    Residue class stores the information associated with
    a particular residue from NMR-STAR loops
    """
    def __init__(self, type):
        self.type = type
        self.atom_chem_shift = {}
        self.atom_chem_shift_err = {}
    # end __init__

    def add_atom(self, atom_type, chem_shift, chem_shift_err=0.2):
        """
        adds atom to a residue.
        Usage: add_atom(atom_type, chem_shift) default chem_shift_err = 0.2
        or   : add_atom(atom_type, chem_shift, chem_shift_err)
        """
        self.atom_chem_shift[atom_type] = chem_shift
        self.atom_chem_shift_err[atom_type] = chem_shift_err
    # end add_atom

    def get_chem_shift(self, atom_type):
        """
        retrieve atom chemical shift from a residue
        """
        if self.atom_chem_shift.get(atom_type):
            return self.atom_chem_shift[atom_type]
        else:
            return None
    # end get_atom

    def get_chem_shift_err(self, atom_type):
        """
        retrieve atom chemical shift error from a residue
        """
        if self.atom_chem_shift_err.get(atom_type):
            return self.atom_chem_shift_err[atom_type]
        else:
            return None
    # end get_atom_err

    def get_resi_type(self):
        return self.type

    def __str__(self):
        s = ''
        form = ('N', 'C', 'CA', 'CB', 'CG', 'CD')
        for item in form:
            if self.atom_chem_shift.get(item):
                s = s + item + ': ' + str(self.atom_chem_shift[item]) + ' '
            else:
                s = s + item + ': ' + '- '
        return s
# end Residue class definition


class Sequence(object):
    """
    Constructs a dictionary of Residues in a sequence
    provides access to them, supports getting slices from 
    sequence and random removal of residues.
    ***
    Residues contain chemical shift information, so 
    Sequence contains only assigned residue objects
    if a particular residue is not assigned it will not 
    be part of the Sequence data structure
    ***
    Instance variables:    
    sequence = {}   - sequence dictionary containing {id : Residue}
    id_list = []    - list of valid ids
    seqAA = []      - amino acid sequence, including ALL residues
                      both assigned/unassigned ones
    """

    def __init__(self, filename = None):
        """
        Constructs a new Sequence object, optionally can be called with 
        srt filename
        """
        self.sequence = {}
        self.id_list = []
        self.seqAA = {}
        self.name = 'bmrb_entry'
        if filename is not None:
            self._init_filename(filename)
    # end __init__
    
    def _init_filename(self, filename):
        """
        parses the srt file using bmrb parser and creates Sequence
        data structures based on NCACX and N(CO)CACX types of experiments
        """
        seq = Sequence()
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
        # it turned out that there could be more than one 
        # assigned_chemical_shifts dataframe - 
        # it is not quite yet clear how I could choose the one i need
        for saveframe in saveframeCS:
            loopCS = saveframe.getLoopByCategory('Atom_chem_shift')
            Data = loopCS.getDataByTag(['Seq_ID','Comp_ID','Atom_ID','Val'])
            for item in Data:
                seq.add_atom(int(item[0]), item[1], item[2], item[3])
            #print seq
            """
            table = assignment.Table(seq)
            table.init_signal_table(table_type = 'N(CO)CACB')    
            table_list.append(table.copy())    
            table.init_signal_table(table_type = 'NCACX')    
            table_list.append(table.copy())
            """            
            #print table
        self.sequence = seq.sequence
        self.id_list = seq.id_list
        self.seqAA = seq.seqAA
        self.name = 'bmrb_entry'
     
    # end read_sequence

    
        
    def change_name(self, name):
        self.name = name
    # change name of the sequence
        
    def __len__(self):
        return len(self.sequence)
    # end
    
    def set_seqAA(self, dictionary):
        """
        sets amino acid sequence, including ALL residues
        both assigned/unassigned ones
        """
        keys = dictionary.keys()
        for id in keys:
            val = dictionary[id]
            if val in _RESI_TYPES_BRIEF:
                self.seqAA[id] = val
            elif val in _RESI_TYPES:
                self.seqAA[id] = _RESI_DICT[val]
            else:
                print "ERROR RESIDUE TYPE"        
    # end set_seqAA
        
    def get_seqAA(self):
        """
        returns amino acid sequence, including ALL residues
        both assigned/unassigned ones
        """        
        return self.seqAA
    # end get_seqAA
        
    def print_seqAA(self):
        """
        returns amino acid sequence as a string of characters
        """
        # return self.sequence.get_seqAA()
        string = ''
        d = self.seqAA
        keys = d.keys()
        keys.sort()
        for i in keys:
            string = string + d[i]
        return string
    # end print_seqAA
        
    def add_atom(self, resn,
                 res_type,
                 atom_type,
                 chem_shift,
                 chem_shift_err=0.2):
        """
        Adds atom to sequence
        resn - residue number
        res_type - type of residue Ala, Met etc.        
        """
        if self.sequence.get(resn) is not None:
            # print "residue exists"
            self.sequence[resn].add_atom(atom_type, chem_shift, chem_shift_err)
        else:
            residue = Residue(res_type)
            residue.add_atom(atom_type, chem_shift, chem_shift_err)
            self.sequence[resn] = residue
    # end add_atom

    def __str__(self):
        s = ''
        keylist = self.sequence.keys()
        # keylist = [int(x) for x in keylist]
        keylist.sort(key=int)
        for key in keylist:
            s = s + "ID: " + str(key)
            s = s + " Residue: " + self.sequence[key].type
            resi = self.sequence[key]
            s = s + " atoms: " + str(resi) + '\n'
        return s
    # end __str__

    def copy(self):
        """
        makes a deep copy of the Sequence object        
        """
        copy = Sequence()
        for item in self.sequence:
            copy.sequence[item] = self.sequence[item]
        copy.seqAA = self.seqAA
        copy.id_list = self.id_list
        return copy
    # end copy

    def get_slice(self, id_list):
        """
        returns a slice from the Sequence based on the list of IDs
        preserves the amino acid sequence
        """
        new = Sequence()
        for i in id_list:
            if self.sequence.get(i):
                new.sequence[i] = self.sequence[i]
                new.id_list.append(i)        
        return new

    def get_id_list(self):
        """
        returns the list of Residue IDs in the Sequence
        """
        if self.id_list == []:
            list = []
            keylist = self.sequence.keys()
            keylist.sort(key=int)
            for key in keylist:
                list.append(key)
            self.id_list = list
            return list
        else:
            return self.id_list

    def delete_id(self, id):
        """
        deletes a Residue with specified ID
        """
        if id in self.id_list:
            del self.sequence[id]
            self.id_list.remove(id) # this works because id is unique
    # end delete_id

    def get_residue(self, id):
        """
        returns a Residue of the Sequence based on its ID
        """
        if id in self.sequence:
            return self.sequence[id]
        else:
            return None
            
    def knock_out(self, threshold=0.1, adj='no'):
        """
        Constructs a new Sequence, where threshold * 100 % of 
        Residue are randomly removed
        Two modes are possible :
        1) a sample of threshold*Sequence_length residues are 
        picked randomly, adjacency (adj)= 'yes'
        2) threshold*Sequence_length residues are picked in succession
        along the sequence adjacency (adj) = 'no'
        """
        import random as rnd
        seq = self
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
# end Sequence class definition
# Table constructs and outputs a table of entries
# corresponding to particular peaks
class Table(object):
    """
    Table class uses Sequence data structure to
    construct tables in MCASSIGN input format
    ***
    Instance variables:
    sequence   -  makes an internal copy of Sequence
    table_type - NCACX, NCOCX etc
    entry_list - list of _Entry class objects containing table entries    
    """
    THIS = True
    THAT = False

    def __init__(self, seq):
        if type(seq).__name__ != 'Sequence':
            sys.exit('No sequence at the input')
        self.sequence = seq.copy()  # ensures immutability of Table
        # when sequence is changed
        self.table_type = ''
        self.entry_list = []
    # end __init__

    def copy(self):
        """
        makes a deep copy of the Table
        """
        copy = Table(self.sequence)
        copy.entry_list = self.entry_list
        copy.table_type = self.table_type
        return copy
    # end copy

    def get_entry(self, id):
        """
        to be implemented later
        """
        pass
    # end get_entry_degeneracy

    def set_entry(self, id, resi_type=[], degeneracy=1):
        """
        to be implemented later
        currently only works to update the
        degeneracy, or residue type

        value could be either integer for degeneracy
        of a string of characters 'YHF' for resi_type

        it currently does not support adding entries that do not 
        yet exist, only changes the entry field
        """
        if type(id) is not int:
            return
            
        if id < 0 or id > len(self.entry_list):
            return
        
        self.entry_list[id].degeneracy = degeneracy
                   
        if resi_type == []:
            return
        else:
            l = []
            for item in resi_type:                
                print item
                if item in _RESI_TYPES:
                    l.append(item)
                if item in _RESI_TYPES_BRIEF:
                    l.append(_RESI_DICT_REVERSE[item])                
                self.entry_list[id].resi_type = l                 
    # end set_entry

    def init_signal_table(self, table_type='NCACX',
                          format=['N', 'C', 'CA', 'CB', 'CG']):
        """
        prepares the MCASSIGN signal table based on the type of table
        and specified output format
                                 
        Currently supported table types
            NCACX
            N(CO)CACB
            NCOCX
        """
        # based on the format split the atoms into two groups
        # ones that belong to one residue and ones that belong to the next
        self.entry_list = []
        this_that = []
        
        self.table_type = table_type
        if table_type == 'NCACX':
            for item in format:                
                this_that.append(self.THIS)
                
        if table_type == 'NCOCX':
            for item in format:
                if item == 'N':
                    this_that.append(self.THAT)
                else:
                    this_that.append(self.THIS)
                    
        if table_type == 'N(CO)CACB':
            for item in format:
                if item == 'N':
                    this_that.append(self.THAT)
                elif item == 'C':
                    this_that.append(None)                    
                else:
                    this_that.append(self.THIS)
        
        # use the prepared lists as inputs ....
        
        format2 = dict(zip(format, this_that))
        id_list = self.sequence.get_id_list()
        for id in id_list:
            resi_this = self.sequence.get_residue(id)
            resi_that = self.sequence.get_residue(id + 1)
            # for non-null residues construct entry
            chem_shift_list = []
            chem_shift_err_list = []
            atom_type_list = []
            if (resi_that and resi_this) is not None:                
                keylist = format2.keys()
                for key in keylist:
                    if format2[key] is self.THIS:
                        cs = resi_this.get_chem_shift(key)
                        cs_err = resi_this.get_chem_shift_err(key)
                    elif format2[key] is self.THAT:
                        cs = resi_that.get_chem_shift(key)
                        cs_err = resi_this.get_chem_shift_err(key)
                    else:
                        cs = None
                        cs_err = None                        
                        pass
                    atom_type_list.append(key)
                    chem_shift_list.append(cs)
                    chem_shift_err_list.append(cs_err)
                # endfor
                
                e = _Entry(atom_type_list,
                           chem_shift_list,
                           chem_shift_err_list,
                           resi_type=[resi_this.get_resi_type()])
                self.entry_list.append(e)
            # endif                    
    # end init_signal_table
                
    def __str__(self):
        string = ''
        count = 0
        for item in self.entry_list:
            string = string + str(item) + "\n"
            count = count + 1          
        string = str(count) + "\n" + string
        return string
        # endfor
        
    # end __str__
        
    def print_seqAA(self):
        """
        returns amino acid sequence as a string of characters
        """
        # return self.sequence.get_seqAA()
        string = ''
        d = self.sequence.get_seqAA()
        keys = d.keys()
        keys.sort()
        for i in keys:
            string = string + d[i]
        return string
    # end print_seqAA
        
    def print_mcassign(self, format=['N', 'CA', 'CB', 'CG']):
        """
        returns a string with a signal table in MCASSIGN format
        """
        output = ""
        # table header (excessive)        
        # output = str(self.table_type) + "\n"
        output = output + \
            str(len(self.entry_list)) + "," +  \
            str(len(format)) + "\n"
        # FROM OLDER VERSION
        # for item in format:
        #    output = output + item + "  "
        # output = output + "\n"
        # table body        
        for entry in self.entry_list:
            chem_shift_list = entry.get_cs(format)
            chem_shift_err_list = entry.get_cs_err(format)
            for item in chem_shift_list:
                if item is None:
                    output = output + "1e6" + "\t"
                else:
                    output = output + str(item) + "\t"
            for item in chem_shift_err_list:
                if item is None:
                    output = output + "1e6" + "\t"
                else:
                    output = output + str(item) + "\t"        
            output = output + str(entry.degeneracy) + "\t"
            
            for item in entry.resi_type:
                output = output + _RESI_DICT[item]
            output = output + "\n"
        return output
        pass
    # end print_MCASSIGN
# end Table class definition
        
class Mcassign(object):
    """
    Data structure constructed from the output tables  of MCASSIGN program.
    """
    def __init__(self, filename):        
        alist = []
        with open(filename, 'r') as f:
            # process each line
            for line in f:
                l = re.split('\s+',line)
                l = [item for item in l if item!='']
                alist.append(l)            
            # end line processing
        self.list_strings = alist
        # end reading        
        
        self.list_floats = []
        # calculate the list of indices with chemical shifts
        length = len(self.list_strings[0])
        s1 = 3
        s2 = 3 + (length - 4)/2 - 1
        s3 = 3 + (length - 4)/2 + 1
        s4 = length - 1
        a = range(s1, s2 + 1, 1)
        print a
        a1 = range(s3, s4 + 1, 1)
        print a1
        a.extend(a1)
        l = []
        for item in self.list_strings:
            item2 = []
            for i in a:
                item2.append(float(item[i]))           
            l.append(item2)        
        self.list_floats = l
        
        self.list_table_line = []
        for item in self.list_floats:
            self.list_table_line.append(_TableLine(item))        
    # end __init__
            
    def compare(self, that):
        """
        Compares two tables and outputs a dictionary
        {'similar' : sim, 'different': diff}
        sim  - table of strings found in table that
        diff - table of strings not found in table that
        """
        similar = []
        different = []
        for line, line2 in zip(self.list_table_line, self.list_strings):
            if line in that.list_table_line:
                similar.append(line2)
            else:
                different.append(line2)
        return {'similar': similar, 'different': different}
        # calculate the list of indices with chemical shifts
    # end compare
# end Mcassign class definition

        


# entries in a Table
class _Entry(object):
    """
    helper class for Table entries
    ***
    Instance variables:
    degeneracy - integer
    resi_type  - residue type
    entry      - dictionary  {atom_type : [chem_shift, chem_shift_err]}    
    """
    
    def __init__(self, atom_type_list,
                 chem_shift_list,
                 chem_shift_err_list=[],
                 resi_type=['None']):
        
        self.degeneracy = 1
        self.resi_type = resi_type
        self.entry = {}
        if chem_shift_err_list == []:
            for itemp in chem_shift_list:
                chem_shift_err_list.append(0.2)
        if len(chem_shift_list) != len(atom_type_list) or \
           len(atom_type_list) != len(chem_shift_err_list):
                sys.exit('Wrong entries: number ' +
                         'of items in the argument lists must be the same')
       
        else:         
            values = []
            values = zip(chem_shift_list, chem_shift_err_list)           
            self.entry = dict(zip(atom_type_list, values))
    # end __init__
    
    def get_cs(self, atom_list):
        """
        return chemical shift for nuclei in the input list
        """
        l = []
        for item in atom_list:
            if item in self.entry:
                l.append(self.entry[item][0])
            else:
                l.append(None)
        return l
    # end get_cs
        
    def get_cs_err(self, atom_list):
        """
        return chemical shift erros for nuclei in the input list
        """
        l = []
        for item in atom_list:
            if item in self.entry:
                l.append(self.entry[item][1])
            else:
                l.append(None)
        return l     
    # end get_cs_err
        
    def __str__(self):
        return str(self.entry)
        # return string
    # end __str__
# end Entry class definition

class _TableLine(object):
    """
    helper class representing the MCASSIGN output table entries
    contains __cmp__ method allowing to compare the entries in 
    different MCASSIGN tables
    """
    def __init__(self, l):
        self.list = l
    # end __init__
        
    def __cmp__(self, that):
        l = []
        for item1, item2 in zip(self.list, that.list):
            l.append(item1 - item2)
        
        # equality testing        
        s = 0
        for item in l:
            s = s + item*item
        if s < 0.1:
            return 0
        
        # determine item inequality
        for item in l:
            if item < 0:
                return -1
            if item > 0:
                return 1        
# end _TableLine helper class to Mcassign

class Run(object):
    def __init__(self, sequence, params):
        self.sequence = sequence
        self.params = params
        table_list = []
        table = Table(sequence)
        table.init_signal_table(table_type = 'N(CO)CACB')    
        table_list.append(table.copy())    
        table.init_signal_table(table_type = 'NCACX')    
        table_list.append(table.copy())
        self.run = {'sequence': self.sequence, 'tables' : 
             {'NCACX' : table_list[1], 'N(CO)CACB' : table_list[0]}}
             
    def start(self):
        sequence_name = self.run['sequence'].name

        self.params['protein_sequence'] = sequence_name + '_seq.txt'
        self.params['connection_table'] = sequence_name + '_conn.txt'
        self.params['signal_tables'] = {'NCACX' : sequence_name + '_NCACX.txt',
                                   'N(CO)CACB' : sequence_name + '_NCOCACB.txt'}
        self.params['connection_table'] = sequence_name + '_conn.txt'
        self.params['output_prefix'] = sequence_name[0:3] + 'dd'
        self.params['output_file'] = sequence_name[0:3] + '.out'
        
        # print "*******SEQUENCE****************"
        # print write_AAsequence(a['sequence'])
        with open(self.params['protein_sequence'], 'w') as f:
            f.write(write_AAsequence(self.run['sequence']))
            file_list = []
            file_list.append(self.params['protein_sequence'])
        
        # print "*******CONNECTIONS*************"
        with open(self.params['connection_table'], 'w') as f:
            f.write(write_connections())
            file_list.append(self.params['connection_table'])
        
        # print "*******NCACX SIGNAL TABLE******"
        with open(self.params['signal_tables']['NCACX'], 'w') as f:
            f.write(self.run['tables']['NCACX'].print_mcassign())
            file_list.append(self.params['signal_tables']['NCACX'])    
        
        #print "*******N(CO)CACB SIGNAL TABLE**"
        with open(self.params['signal_tables']['N(CO)CACB'], 'w') as f:
            f.write(self.run['tables']['N(CO)CACB'].print_mcassign())
            file_list.append(self.params['signal_tables']['N(CO)CACB'])    
        
        with open(sequence_name + '.in', 'w') as f:
            f.write(write_input(self.params))
            input_file = sequence_name + '.in'
            file_list.append(input_file)
            
        file_list.append(self.params['output_file'])
        log_file = sequence_name + '.log'
        file_list.append(log_file)
        
        # RUNNING MCASSIGN WITH INPUT FILES PREPARED ABOVE
        print "Starting..."
        print "mcassign2b.exe < " + input_file + " > " + log_file
        sys.stdout.flush()
        os.system("mcassign2b.exe < " + input_file + " > " + log_file )
            
        #   BOOKKEEPING AFTER CALCULATION IS COMPLETE
        import glob
        file_list.extend(glob.glob(self.params['output_prefix'] + '*'))    
        
        try:
            os.makedirs(self.run['sequence'].name)
        except OSError:
            pass
        
        import shutil
        for file in file_list:
            try:
                shutil.move(file, self.run['sequence'].name + '/' + file)
            except IOError:
                print file + " does not exist"
                            
# end Run class

# various helper procedures
def write_AAsequence(sequence):
    """
    returns a string with amino acid characters
    """
    return sequence.print_seqAA() + '\n'
    
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

# test client
"""

# TESTING RESIDUE CLASS
a = Residue('Ala')
a.add_atom('CA', 44.5)
form = ['CA', 'CB', 'N']
print "this is a formatted output:\n"
print a
print a.get_chem_shift('CA')
if a.get_chem_shift('N') is None:
    print "no such item"

# TESTING SEQUENCE CLASS
print "*****testing Sequence****"
b = Sequence()
b.add_atom(1, 'Cys', 'CA', 34.5)
b.add_atom(2, 'Ala', 'CA', 5.5)
b.add_atom(2, 'Ala', 'CB', 6.5)
b.add_atom(2, 'Ala', 'CG', 7.5)

b.add_atom(3, 'Bla', 'CA', 5.5)
b.add_atom(3, 'Bla', 'CB', 6.5)
b.add_atom(3, 'Bla', 'CG', 7.5)

c = b.copy()
d = c.get_slice([1, 2, 3])
print d

# print b.sequence
residue = b.get_residue(4)
if residue is None:
    print "None"
else:
    print residue

print b.get_id_list()

# TESTING _ENTRY CLASS
print "****testing _Entry*****"
entry = _Entry(['CA', 'CB', 'CD'], [1, 2, 3])
print entry
#entry = _Entry(['CA', 'CB'], [1, 2, 3])
print entry
print "+++++++++++++++++"
print entry.get_cs(['CA', 'CB', 'CG'])
print "+++++++++++++++++"


# TESTING TABLE CLASS
print "*******testing Table*********"

atom_names = ('CA', 'CB', 'CG', 'CD', 'CO', 'N')
resi_types = ('Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly',
              'His', 'Ile', 'Lys', 'Leu', 'Met', 'Gln')

import random

random.seed(a = 1)

seq1 = Sequence()
for i in range(10):    
    resi_name = random.choice(resi_types)
    for name in atom_names:
        seq1.add_atom(i + 1, resi_name, name, round((random.gauss(50, 20)), 2))

print "**********PRINT SEQUENCE**********"
print seq1
table2 = Table(seq1)
table2.init_signal_table()
print "**********PRINT NCACX TABLE**********"
print table2
table2.init_signal_table(table_type = 'N(CO)CACB')
print "**********PRINT N(CO)CACB TABLE**********"
print table2
"""
