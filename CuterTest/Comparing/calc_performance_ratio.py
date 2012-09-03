#!/usr/bin/python

import subprocess
import sys
import re

if len(sys.argv) < 2:
  print "Precisa do prefixo"
  exit()

file_dci = open(sys.argv[1] + '.dcicpp','r')
file_ipopt = open(sys.argv[1] + '.ipopt','r')
file_out = open(sys.argv[1] + '.ratio','w')

maxvalue = 'max'

for line_dci in file_dci:
  line_ipopt = file_ipopt.readline()

  problem = line_dci.split(' ')[0]

  result_dci = line_dci.split(' ')[1]
  result_ipopt = line_ipopt.split(' ')[1]
  line_dci = float(line_dci.split(' ')[2])
  line_ipopt = float(line_ipopt.split(' ')[2])

  m = min(line_dci, line_ipopt)
  if result_dci != 'Converged':
    line_dci = maxvalue
  elif m == 0 and line_dci > 0:
    line_dci = maxvalue
  elif m == 0 and line_dci == 0:
    line_dci = 1
  else:
    line_dci = line_dci/m

  if result_ipopt != 'Optimal':
    line_ipopt = maxvalue
  elif m == 0 and line_ipopt > 0:
    line_ipopt = maxvalue
  elif m == 0 and line_ipopt == 0:
    line_ipopt = 1
  else:
    line_ipopt = line_ipopt/m

  file_out.write(' '.join(str(x) for x in [problem, line_dci, line_ipopt,'\n']))

file_dci.close()
file_ipopt.close()
