#!/usr/bin/python

# Two files are given as arguments. We must have both files
# with the same problems in the same order.

import subprocess
import sys
import re
import math

if len(sys.argv) < 3:
  print "Need two latex_* files for comparing"
  exit()

file1 = open(sys.argv[1], 'r');
file2 = open(sys.argv[2], 'r');

bf1 = 0;
bf2 = 0;
bt1 = 0;
bt2 = 0;

print "Prob. Name %objfun% %f_norm% %%time%% %time_norm%"

for line1 in file1:
  line1 = line1.split();
  line2 = file2.readline();
  line2 = line2.split();

  if (line1[0] != line2[0]):
    print "ERROR: Files must have the same problems in the same order"
    exit()
  
  f1 = float(line1[6]);
  f2 = float(line2[6]);
  t1 = float(line1[14]);
  t2 = float(line2[14]);

  if (f1 < f2):
    bf1 += 1
  elif (f2 < f1):
    bf2 += 1

  if (t1 < t2):
    bt1 += 1
  elif (t2 < t1):
    bt2 += 1

  f = f1 - f2;
  fnorm = f/(min(abs(f1),abs(f2)) + 1)
  t = t1 - t2
  tnorm = t/(min(abs(t1),abs(t2)) + 1)

  print ' '.join(str(x) for x in 
          [line1[0].ljust(10), "%10.3e" % f, 
            "%10.3e" % fnorm, 
            "%10.3e" % t, 
            "%10.3e" % tnorm])


print ' '.join(str(x) for x in ['best f1 =', bf1])
print ' '.join(str(x) for x in ['best f2 =', bf2])
print ' '.join(str(x) for x in ['best t1 =', bt1])
print ' '.join(str(x) for x in ['best t2 =', bt2])
