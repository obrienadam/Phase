#!/usr/bin/env python2

import glob
import os, sys
import tarfile

for fname in (fname for arg in sys.argv[1:] for fname in glob.glob(arg)):
  with tarfile.open(fname, 'r') as f:
    print 'Extracting archive "{}"...'.format(fname)
    f.extractall(path=os.path.dirname(fname))
