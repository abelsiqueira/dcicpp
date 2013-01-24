#!/usr/bin/python
# IN: cond.something. 
# OUT: something.wellcond.list
#
# Removes entries from IN ill conditioned

import subprocess
import sys

ill_max = 1e+4

if len(sys.argv) < 2:
  print "Need argument"
  exit()
  
f = open(sys.argv[1], 'r')

for line in f:
  line = line.split(' ');
  len_line = len(line);
  try:
    ill_value = float(line[len(line)-1]);
  except ValueError:
    continue
  if ill_value < ill_max:
    problem = line[1];
    problem = problem[0:problem.find('_')];
    print problem;
