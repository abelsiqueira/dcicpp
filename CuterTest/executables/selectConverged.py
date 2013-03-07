#!/usr/bin/python

import subprocess
import sys
import re

if len(sys.argv) < 3:
  print "Precisa do prefixo dos testes e da saida"
  print "Usage: ./selectConverged testIn testOut [strict]"
  exit()

if len(sys.argv) < 4:
  strict = True;
else:
  if sys.argv[3] == 'True':
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
file_out_dci = open(sys.argv[2] + '.dcicpp','w')
file_out_ipopt = open(sys.argv[2] + '.ipopt','w')

for line_dci in file_dci:
  line_ipopt = file_ipopt.readline()

  st_dci = line_dci.split(' ')[1]
  st_ipopt = line_ipopt.split(' ')[1]

  if st_dci in dci_results and st_ipopt in ipopt_results:
    file_out_dci.write(line_dci)
    file_out_ipopt.write(line_ipopt)

file_dci.close()
file_ipopt.close()
file_out_dci.close()
file_out_ipopt.close()
