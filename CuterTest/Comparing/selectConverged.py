#!/usr/bin/python

import subprocess
import sys
import re

if len(sys.argv) < 3:
  print "Precisa do prefixo dos testes e da saida"
  print "Usage: ./selectConverged testIn testOut"
  exit()

file_dci = open(sys.argv[1] + '.dcicpp','r')
file_ipopt = open(sys.argv[1] + '.ipopt','r')
file_out_dci = open(sys.argv[2] + '.dcicpp','w')
file_out_ipopt = open(sys.argv[2] + '.ipopt','w')

for line_dci in file_dci:
  line_ipopt = file_ipopt.readline()

  status = line_dci.split(' ')[1]

  if status == 'Converged':
    file_out_dci.write(line_dci)
    file_out_ipopt.write(line_ipopt)

file_dci.close()
file_ipopt.close()
file_out_dci.close()
file_out_ipopt.close()
