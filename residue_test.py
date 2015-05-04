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
"""
import sys


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

    def add_atom(self, atom_type, chem_shift, chem_shift_err = 0.2):
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

    def __str__(self):
        s = ''
        form = ('N', 'CO', 'CA', 'CB', 'CG', 'CD')
        for item in form:
            if self.atom_chem_shift.get(item):
                s = s + item + ': ' + str(self.atom_chem_shift[item]) + ' '
            else:
                s = s + item + ': ' + '- '
        return s
# end Residue class definition

"""
Sequence contains residues in a dictionary
and provides access to the elements of sequence
"""
class Sequence(object):
    """
    sequence = {}
    id_list = []
    """
    
    def __init__(self):
        self.sequence = {}
        self.id_list = []
    # end __init__

    def add_atom(self, resn,
                 res_type,
                 atom_type,
                 chem_shift,
                 chem_shift_err = 0.2):
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

    def get_residue(self, id):
        if id in self.sequence:
            return self.sequence[id]
        else:
            return None
# end Sequence class definition


# Table constructs and outputs a table of entries
# corresponding to particular peaks
class Table(object):
    """
    entry_list = []
    sequence = Sequence()
    """
    THIS = True;
    THAT = False;

    def __init__(self, seq):
        if type(seq).__name__ != 'Sequence':
            sys.exit('No sequence at the input')
        self.sequence = seq
        self.entry_list = []
    # end __init__

    def set_entry_degeneracy(self,
                             name_list,
                             cs_list,
                             id,
                             degeneracy,
                             type_list):
        pass

    def set_entry_type(self, id, type):
        pass

    def init_signal_table(self, table_type = 'NCACX',
                                 format = ['N', 'CO', 'CA', 'CB', 'CG']):
        """
        this method prepares the signal table based on the type of table
        and specified output format
        Usage: print_signal_table(table_type = 'NCACX',
                                 format= ['N', 'CO', 'CA', 'CB', CG])
        """
        # based on the format split the atoms into two groups
        # ones that belong to one residue and ones that belong to the next
        self.entry_list = []
        this_that = []
        
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
                elif item == 'CO':
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
                e = _Entry(atom_type_list, chem_shift_list, chem_shift_err_list)
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
# end Table class definition


# entries in a Table
class _Entry(object):
    # self:
    # entry = {}
    def __init__(self, atom_type_list,
                       chem_shift_list,
                       chem_shift_err_list = []):
        
        self.entry = {}
        if chem_shift_err_list == []:
            for itemp in chem_shift_list:
                chem_shift_err_list.append(0.2)
        if len(chem_shift_list) != len(atom_type_list) or \
             len(atom_type_list) != len(chem_shift_err_list):
                sys.exit('Wrong entries: number ' +
                         'of items in the argument lists must be the same')
       
        else:
            self.entry = dict(zip(atom_type_list, chem_shift_list))
    # end __init__

    def get(self, name_list):
        return self.entry[name_list]
    

    def __str__(self):
        return str(self.entry)
        #return string
# end Entry class definition


# test client
a = Residue('Ala')
a.add_atom('CA', 44.5)
form = ['CA', 'CB', 'N']
print "this is a formatted output:\n"
print a
print a.get_chem_shift('CA')
if a.get_chem_shift('N') is None:
    print "no such item"

# test Sequence class methods
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

# test  _Entry class
print "****testing _Entry*****"
entry = _Entry(['CA', 'CB', 'CD'], [1, 2, 3])
print entry
#entry = _Entry(['CA', 'CB'], [1, 2, 3])
print entry

print "*******testing Table*********"

table = Table(b)
table.init_signal_table()
print table


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

print seq1
table2 = Table(seq1)
table2.init_signal_table()
print table2
table2.init_signal_table(table_type = 'N(CO)CACB')
print table2

