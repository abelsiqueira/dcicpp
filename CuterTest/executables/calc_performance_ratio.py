#!/usr/bin/python

import subprocess
import sys
import re

if len(sys.argv) < 2:
  print "Precisa do prefixo"
  exit()

if len(sys.argv) < 3:
  strict = True;
else:
  if sys.argv[2] == 'True':
    strict = True
  else:
    strict = False

if not strict:
  dci_results = ['Converged','stationary']
  ipopt_results = ['Optimal','Converged','Solved']
else:
  dci_results = ['Converged']
  ipopt_results = ['Optimal','Solved']

file_dci = open(sys.argv[1] + '.dcicpp','r')
file_ipopt = open(sys.argv[1] + '.ipopt','r')
file_out = open(sys.argv[1] + '.ratio','w')

smallValue = 1e-2
maxvalue = 'max'

for line_dci in file_dci:
  line_ipopt = file_ipopt.readline()

  problem = line_dci.split(' ')[0]
  if (problem != line_ipopt.split(' ')[0]):
    print "NAMES DONT MATCH"
    exit(1)

  result_dci = line_dci.split(' ')[1]
  result_ipopt = line_ipopt.split(' ')[1]
  line_dci = line_dci.split(' ')[2]
  if line_dci != '':
    line_dci = float(line_dci)
  else:
    line_dci = 0
  line_ipopt = line_ipopt.split(' ')[2]
  if line_ipopt.strip() != '':
    line_ipopt = float(line_ipopt)
  else:
    line_ipopt = 2

  m = min(line_dci, line_ipopt)
  if result_dci not in dci_results:
    line_dci = maxvalue
  elif m == 0 and line_dci > 0:
    line_dci = line_dci/smallValue
  elif m == 0 and line_dci == 0:
    line_dci = 1
  else:
    line_dci = line_dci/m

  if result_ipopt not in ipopt_results:
    line_ipopt = maxvalue
  elif m == 0 and line_ipopt > 0:
    line_ipopt = line_ipopt/smallValue
  elif m == 0 and line_ipopt == 0:
    line_ipopt = 1
  else:
    line_ipopt = line_ipopt/m

  file_out.write(' '.join(str(x) for x in [problem, line_dci, line_ipopt,'\n']))

file_dci.close()
file_ipopt.close()
