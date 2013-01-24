#!/bin/bash
# Given $1 and $2, finds
# all files in $1 which are not in $2 (using grep)
# Useful when some test was not made, and you need
# to find it.

for problem in $(cat $1)
do
  grep_output="$(grep $problem *)"
  if [ -z "$grep_output" ];
  then
    echo "$problem not found"
  fi
done
