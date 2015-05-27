# -*- coding: utf-8 -*-
"""
Created on Mon May 25 17:12:43 2015
sfafsf
@author: Alexey
MCASSIGN program parameters
"""
params = {}
params['protein_sequence'] = \
                           ''
params['number_of_signal_tables'] = \
                           2
params['signal_tables'] = \
                           {'NCACX': 'NCACX.txt', 'N(CO)CACB' :'NCOCX.txt'}
params['connection_table'] = \
                           ''
params['random_seed'] = \
                           '1' # if 0 the seed will derived from system clock
params['good_connections'] = \
                           {'initial' : 0, 'final': 10}
params['bad_connections'] = \
                           {'initial' : 0, 'final': 20}
params['edge_penalty'] = \
                           {'initial' : 0, 'final': 4}
params['signal_usage'] = \
                           {'initial' : 0, 'final': 1}
params['step_number'] = \
                           {'weighting' : 20, 'Monte Carlo': 100000}
params['runs_number'] = \
                           10
params['output_prefix'] = \
                           ''  # 5 characters limit
params['output_file'] = \
                           ''
params['initial_assignments'] = \
                            0 # 1 - specified, 0 - not specified
