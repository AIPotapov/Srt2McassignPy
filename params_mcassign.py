# -*- coding: utf-8 -*-
"""
Created on Mon May 25 17:12:43 2015
@author: Alexey
MCASSIGN program parameters
it is just a better layout of what comes into MCASSIGN program as parameters
typed from the command line
"""
default = {}
default['protein_sequence'] = \
                           ''
default['number_of_signal_tables'] = \
                           2
default['signal_tables'] = \
                           {'NCACX': 'NCACX.txt', 'N(CO)CACB' :'NCOCX.txt'}
default['connection_table'] = \
                           ''
default['random_seed'] = \
                           '0' # if 0 the seed will derived from system clock
default['good_connections'] = \
                           {'initial' : 0, 'final': 10}
default['bad_connections'] = \
                           {'initial' : 0, 'final': 20}
default['edge_penalty'] = \
                           {'initial' : 0, 'final': 4}
default['signal_usage'] = \
                           {'initial' : 0, 'final': 1}
default['step_number'] = \
                           {'weighting' : 20, 'Monte Carlo': 100000}
default['runs_number'] = \
                           10
default['output_prefix'] = \
                           ''  # 5 characters limit
default['output_file'] = \
                           ''
default['initial_assignments'] = \
                            0 # 1 - specified, 0 - not specified
