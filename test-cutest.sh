#!/bin/bash

set -v

tmpdir=tmpdir
mkdir -p $tmpdir

cp Tests/dcicpp.spc $tmpdir
cd $tmpdir
c=0
for n in $(seq 1 119)
do
  g=$(rundcicpp -D HS$n -lgfortran -lgfortranbegin | grep Converged)
  if [ ! -z "$g" ]; then
    c=$(($c+1))
  fi
done

echo "Convergence count: $c"
[ $c -lt 110 ] && exit 1
