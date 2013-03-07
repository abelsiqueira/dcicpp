#!/usr/bin/python

import subprocess
import sys
import re

if len(sys.argv) < 2:
  print "arg"
  exit()

p = subprocess.Popen('wc -l ' + sys.argv[1] + '.ratio',
    stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = p.communicate()
out = out.split()
np = int(out[0])

t_list = open(sys.argv[1] + '.t_list','r')

print ''.join(str(x) for x in ['t = zeros(',np,',1);'])
print ''.join(str(x) for x in ['dci_f = zeros(',np,',1);'])
print ''.join(str(x) for x in ['ipopt_f = zeros(',np,',1);'])

k = 1
for t in t_list:
  t = float(t)
  print ''.join(str(x) for x in ['t(',k,") = ", t,';'])

  dci_f = 0.0
  ipopt_f = 0.0
  file = open(sys.argv[1] + '.ratio','r')
  for line in file:
    line = line.split()
    if float(line[1]) <= t:
      dci_f += 1
    if float(line[2]) <= t:
      ipopt_f += 1
  file.close()
  dci_f = float(dci_f/float(np))
  ipopt_f = float(ipopt_f/np)
  print ''.join(str(x) for x in ['dci_f(',k,") = ", dci_f,';'])
  print ''.join(str(x) for x in ['ipopt_f(',k,") = ", ipopt_f,';'])
  k = k + 1
t_list.close()

p = subprocess.Popen('./find_max.py ' + sys.argv[1] + '.t_list',
    stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = p.communicate()
out = float(out)

print 't = sort(t);'
print 'dci_f = sort(dci_f);'
print 'ipopt_f = sort(ipopt_f);'
print 'hold off'
print "semilogx(t,dci_f,'-xb')"
print 'hold on'
print "semilogx(t,ipopt_f,'-or')"
print "legend('dcicpp','ipopt','location','southeast')"
print ''.join(str(x) for x in ["xl = xlim; xlim([1, ", out, "]);"])
print "set(gca, 'YTick', [0:0.1:1])"
print "ylim([0,1]);"
print ''.join(str(x) for x in ["print('",sys.argv[1],".png','-color')"])
