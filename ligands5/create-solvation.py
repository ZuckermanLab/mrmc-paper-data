#!/usr/bin/env python

import sys
#to install parmed: python setup.py install --home=~
#To use parmed: echo $PYTHONPATH=~/lib/python:$PYTHONPATH
from parmed.amber import *
#import parmed
parm = LoadParm(sys.argv[1])
typestart = int(sys.argv[2])
classstart = int(sys.argv[3])



