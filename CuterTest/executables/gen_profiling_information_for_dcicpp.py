#!/usr/bin/python

import subprocess
import sys
import re

def run_cuter (prob):
  cmd = ["runcppcuter -p dcicpp -D " , prob , " > garbage"]
  cmd = ''.join(cmd);
  p = subprocess.Popen(cmd,
      stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  out, err = p.communicate()
  cmd = "grep EXIT garbage; grep 'Setup time' garbage; grep 'Execution time' garbage",
  q = subprocess.Popen(cmd,
      stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  out, err = q.communicate()

  return out

if len(sys.argv) < 2:
  print "Need list"
  exit()

f = open(sys.argv[1],'r');
for problem in f:
  problem = problem.rstrip()
  out = run_cuter(problem)
  lines = out.split('\n');
  result = lines[0].split(' ')[4]
  t1 = float(re.sub('\s+', ' ', lines[1]).split(' ')[2])
  t2 = float(re.sub('\s+', ' ', lines[2]).split(' ')[2])

  print problem, result, t1+t2
