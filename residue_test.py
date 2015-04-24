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
    atomCS = {}  # a dictionary of key:value where key -  atom type
    # such as Ca, Cb, Cg etc and value - chemical shift value
    
    def __init__(self, type):
        self.type = type        
        self.atomCS = {}
    
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
        

class Sequence:
    sequence = {}
    def __init__(self):
        self.sequence = {}
            
    def addAtom(self, resN, resType, atomType, chemShift):
        if self.sequence.get(resN) != None:
            # print "residue exists"
            self.sequence[resN].addAtom(atomType, chemShift)
        else:
            residue = Residue(resType)
            residue.addAtom(atomType, chemShift)
            self.sequence[resN] = residue
            
    def toString(self, format):
        s = ''
        for key in self.sequence:
            s =  s + "ID: " + str(key)
            s = s +  " Residue: " + self.sequence[key].type
            resi = self.sequence[key]
            s =  s + " atoms: " + resi.toString(format) +'\n'
        return s
        
    def copy(self):
        copy = Sequence()
        for item in self.sequence:
            copy.sequence[item] = self.sequence[item]
        return copy
        
    def getSlice(self, list):
        new = Sequence()        
        for i in list:
            if self.sequence.get(i):
                new.sequence[i] = self.sequence[i]
        return new
        
# test client
"""

a = Residue('Ala')
a.addAtom('CA', 44.5)
form = ['CA', 'CB', 'N']
print "this is a formatted output:\n" + a.toString(form)
print a.getAtom('CA')
if a.getAtom('N') == None:
    print "no such item"
    
    
print "*****testing Sequence****"
b = Sequence()
b.addAtom(1, 'Cys', 'CA', 34.5)
b.addAtom(2, 'Ala', 'CA', 5.5)
b.addAtom(2, 'Ala', 'CB', 6.5)
b.addAtom(2, 'Ala', 'CG', 7.5)

b.addAtom(3, 'Bla', 'CA', 5.5)
b.addAtom(3, 'Bla', 'CB', 6.5)
b.addAtom(3, 'Bla', 'CG', 7.5)



c = b.copy()
d = c.getSlice([1,2,3])
format = ('CB','CA','CG','CD')
s =  d.toString(format)
print s
#print b.sequence
"""
