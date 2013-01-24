#!/usr/bin/python
# IN: a file from which you want the first column
# OUT: a list with the first column
#

import subprocess
import sys

if len(sys.argv) < 2:
  print "Need file"
  exit()
  
f = open(sys.argv[1], 'r')

for line in f:
  line = line.strip().split(' ');
  print line[0]
