#!/usr/bin/python

import subprocess
import sys
import re
import math

if len(sys.argv) < 2:
  print "Need prefix for wc files"
  exit()

prefix = sys.argv[1];

cmd = ''.join(["ls ", prefix, "*wc"])
p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, shell=True)
out, err = p.communicate()

out = out.split('\n')

conv    = 0;
maxiter = 1;
maxrest = 2;
rhomax  = 3;
time    = 4;
infeas  = 5;
unlim   = 6;
total   = 7;
exits = ['latex_convergence',
         'latex_maxiter',
         'latex_maxrest',
         'latex_rhomax',
         'latex_timelimit',
         'latex_infeasible',
         'latex_unlimited'
         ];
titles = ["Convergiu",
          "M\\'aximo de Itera\\c{c}\\~oes",
          "M\\'aximo de Restaura\\c{c}\\~oes",
          "$\\rhomax$ pequeno",
          "Limite de Tempo",
          "Estacion\\'ario da Inviabilidade",
          "Ilimitado",
          "Total"
          ];

unc    = 0;
equ    = 1;
ineq   = 2;
gencon = 3;
total  = 4;

types  = ['unc',
          'equ',
          'ineq',
          'gencon'
          ];

table = [[0 for x in xrange(8)] for x in xrange(5)]

for test in out:
  if test == '':
    continue

  for i in range(4):
    if types[i] in test:
      thistype = i;
      break
  else:
    print "ERROR", test
    exit()
  
  f = open(test, 'r');
  for line in f:
    line = line.strip().split()
    for i in range(7):
      if exits[i] in line[1]:
        table[thistype][i] += int(line[0]);
        break;

for i in range(7):
  table[4][i] = 0;
  for j in range(4):
    table[4][i] += table[j][i];

for i in range(5):
  table[i][7] = sum(table[i][0:7])

# Now the printing
print "\\begin{table}[H]"
print "\\centering"
print "\\scriptsize"
print "\\begin{tabular}{|c||c|c||c|c||c|c||c|c||c|c|} \\hline"
print "\\multirow{2}{*}{\\bf Sa\\'ida do Algoritmo} &"
print "\\multicolumn{2}{|c||}{\\bf Irrestritos} &"
print "\\multicolumn{2}{|c||}{\\bf Igualdade} &"
print "\\multicolumn{2}{|c||}{\\bf Desigualdade} &"
print "\\multicolumn{2}{|c||}{\\bf Restri\\c{c}\\~oes Gerais} &"
print "\\multicolumn{2}{|c|}{\\bf Total} \\\\ \\cline{2-11}"
print "& {\\bf \\No} & {\\bf \\%}"
print "& {\\bf \\No} & {\\bf \\%}"
print "& {\\bf \\No} & {\\bf \\%}"
print "& {\\bf \\No} & {\\bf \\%}"
print "& {\\bf \\No} & {\\bf \\%}"
print "\\\\ \\hline"

for ex in range(8):
  print "{\\bf ", titles[ex], " } "
  for ty in range(5):
    d = int(table[ty][ex])
    v = float(d)
    T = float(table[ty][7])
    print "&", "%3d" % d, " &", "%6.2f" % (100*v/T), " ",
  print '\\\\ \\hline'

  
print "\\end{tabular}"
print "\\caption{ caption }"
print "\\label{tab:label}"
print "\\end{table}"
