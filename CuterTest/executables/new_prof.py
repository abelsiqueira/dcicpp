#!/usr/bin/python

import subprocess
import sys
import re

prefix = "clean.out";
name_str = prefix + ".name";
exit_str = prefix + ".exit";
setup_str = prefix + ".setup.time";
exec_str = prefix + ".exec.time";

name_file = open(name_str, 'r');
exit_file = open(exit_str, 'r');
setup_file = open(setup_str, 'r');
exec_file = open(exec_str, 'r');

for problem in name_file:
  exitflag = exit_file.readline();
  setuptime = setup_file.readline();
  exectime = exec_file.readline();

  problem = problem.strip().split()[1];
  exitflag = exitflag.split()[4];
  setuptime = re.sub('\s+', ' ', setuptime).split()[2];
  exectime = re.sub('\s+', ' ', exectime).split()[2];
  time = float(setuptime) + float(exectime);

  print problem, exitflag, time
