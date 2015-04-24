# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 22:26:14 2015

@author: Alexey
Definition and test client for Residue class

other possible modifications - perhaps in future there is going
to be a need for fixed format of atom and residue names
"""



# Residue class stores the information associated with a particular
# residue from NMR-STAR loops
class Residue:
    type = ''  # residue type such as Ala, Cys, Leu etc.
    id = 0     # residue number in the sequence
    atomCS = {}  # a dictionary of key:value where key -  atom type
    # such as Ca, Cb, Cg etc and value - chemical shift value
    
    def __init__(self, type, id):
        self.type = type
        self.id = id
    
    def addAtom(self, atomType, chemShift):
        self.atomCS[atomType] = chemShift      

    def getAtom(self, atomType):
        if self.atomCS.get(atomType):
            return self.atomCS[atomType]
        else:
            return None            
        
    def toString(self, form):
        s = ''
        for item in form:
            if (self.atomCS.get(item)):
                s = s + item + ': ' + str(self.atomCS[item]) + ' '
            else: 
                s = s + item + ': ' + '- '
        return s
        
        
# test client

a = Residue('Ala', 4)
a.addAtom('CA', 44.5)
form = ['CA', 'CB', 'N']
print "this is a formatted output:\n" + a.toString(form)
print a.getAtom('CA')
if a.getAtom('N') == None:
    print "no such item"
    
#print a.toString()        
        
   
