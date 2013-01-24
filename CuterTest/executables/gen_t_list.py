#!/usr/bin/python

import subprocess
import sys

if len(sys.argv) < 2:
  print "ERROR"
  exit()

file = open(sys.argv[1],'r')

for line in file:
  line = line.split();
  print line[1]
  print line[2]
