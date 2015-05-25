# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 22:26:14 2015

@author: Alexey Potapov
Definition and test client for Residue class

other possible modifications - perhaps in future there is going
to be a need for fixed format of atom and residue names

some methods are not safe because no format or null
checking takes place

3 may added get_residue
and get_id_list methods
implementation of Table Class - not tested yet

4 modify Residue class to include chem_shift_err
first working version of the Table class
it not yet thoroughly tested and outline may change in the future

need to work on documentation to the classes

14 may
update the implementation of _Entry class

25 may
move Mcassign and its helper class _TableLine from srt2mcassign

"""
import sys
import re

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
    Residue class stores the information associated with a particular residue
    from NMR-STAR loops
    """
    """
    type = ''  # residue type such as Ala, Cys, Leu etc.
    atom_chem_shift = {}  # a dictionary of key:value where key -  atom type
    atom_chem_shift_err = {}
    """
    # such as Ca, Cb, Cg etc and value - chemical shift value

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
    Sequence contains residues in a dictionary
    and provides access to the elements of sequence

    Instance variables:    
    sequence = {}
    id_list = []
    seqAA = []
    """

    def __init__(self):
        self.sequence = {}
        self.id_list = []
        self.seqAA = {}
    # end __init__
        
    def __len__(self):
        return len(self.sequence)
    # end
    
    def set_seqAA(self, dictionary):        
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
        return self.seqAA
    # end get_seqAA

    def add_atom(self, resn,
                 res_type,
                 atom_type,
                 chem_shift,
                 chem_shift_err=0.2):
        """
        Adds atom to sequence
        Usage: add_atom(self, resn,
                 res_type,
                 atom_type,
                 chem_shift,
                 chem_shift_err = 0.2)
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
        copy = Sequence()
        for item in self.sequence:
            copy.sequence[item] = self.sequence[item]
        copy.seqAA = self.seqAA
        copy.id_list = self.id_list
        return copy
    # end copy

    def get_slice(self, list):
        new = Sequence()
        for i in list:
            if self.sequence.get(i):
                new.sequence[i] = self.sequence[i]
        return new

    def get_id_list(self):
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
        if id in self.id_list:
            del self.sequence[id]
    # end delete_id

    def get_residue(self, id):
        if id in self.sequence:
            return self.sequence[id]
        else:
            return None
            
            
            
            
    def knock_out(self, threshold=0.1, adj='no'):
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
    entry_list = []
    sequence = Sequence()
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
        prepares the signal table based on the type of table
        and specified output format
        Usage: print_signal_table(table_type = 'NCACX',
                                 format= ['N', 'CO', 'CA', 'CB', CG])
                                 
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
        prints the table in MCASSIGN format
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


# entries in a Table
class _Entry(object):
    # self:
    # entry = {}
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
        l = []
        for item in atom_list:
            if item in self.entry:
                l.append(self.entry[item][0])
            else:
                l.append(None)
        return l
    # end get_cs
        
    def get_cs_err(self, atom_list):
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

class Mcassign:
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

# help class for Mcassign class
class _TableLine:
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
# end Mcassing class definition





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
