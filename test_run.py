# -*- coding: utf-8 -*-
"""
Created on Thu May 28 19:55:55 2015

@author: Alexey
"""
import os
import sys
sys.path.append("../sans_python/")
import bmrb
sys.path.append("./")
import assignment
from params_mcassign import *

S = assignment.Sequence('15156.srt')
S.change_name('test')
m = assignment.Run(S, default)
m.start()