#!/usr/bin/env python2

import sys
import os

try:
    directory = sys.argv[1]
except IndexError:
    print "Error: must specify the case name"
    sys.exit(-1)

if not os.path.exists(directory):
    os.makedirs(directory)
else:
    print "directory already exists"
    sys.exit(-1)
    
print 'Creating new case "{}"...'.format(directory)
    
os.chdir(directory)
os.makedirs('case')
os.makedirs('solution')

os.chdir('case')

with open('case.info', 'w') as f:
    f.write('; Case input file for {}\n\n'.format(directory))
    f.write('CaseName {}\n\n'.format(directory))
    f.write('Solver\n{\n\n}\n\n')
    f.write('Properties\n{\n\n}\n\n')
    f.write('Grid\n{\n\n}\n\n')
    f.write('Viewer\n{\n\n}\n\n')
    
with open('boundaries.info', 'w') as f:
    f.write('; Boundary condition input file for {}\n\n'.format(directory))
    f.write('Boundaries\n{\n\n}\n\n')
    
with open('initialConditions.info', 'w') as f:
    f.write('; Initial condition input file for {}\n\n'.format(directory))
    f.write('InitialConditions\n{\n\n}\n\n')
    
print 'Finished setting up case "{}"'.format(directory)