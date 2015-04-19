# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 22:31:26 2015

@author: Alexey
Using bmrb.py to read and process NMR-STAR files
"""


"""
Testing the bmrb module
"""


import sys
sys.path.append("../sans_python/")
import bmrb

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


