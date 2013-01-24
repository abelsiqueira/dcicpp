#!/usr/bin/python
import sys

if len(sys.argv) < 2:
  print "Need arg"
  exit()

file = open(sys.argv[1])

Max = -10000
for line in file:
  thismax = max(float(x) for x in line.split(' '))
  if thismax > Max:
    Max = thismax
print Max
