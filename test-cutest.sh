#!/bin/bash

set -v

tmpdir=tmpdir
rm -rf $tmpdir
mkdir -p $tmpdir

cp tests/dcicpp.spc $tmpdir
cp tests/fast.list $tmpdir
cd $tmpdir
rm -f fail.list
c=0
T=0
for problem in $(cat fast.list)
do
  echo "Running problem $problem"
  g=$(rundcicpp -D $problem -lgfortran -lgfortranbegin | grep Converged)
  if [ ! -z "$g" ]; then
    c=$(($c+1))
  else
    echo $problem >> fail.list
  fi
  T=$(($T+1))
  echo "Partial count: $c/$T = $(echo "scale=2;100*$c/$T"|bc)%"
done

p=$(echo "scale=2;100*$c/$T"|bc)
echo "Convergence: $c/$T = $p%"
if [ $(echo "$p > 90" | bc) -eq 1 ];
then
  exit 0
else
  exit 1
fi
