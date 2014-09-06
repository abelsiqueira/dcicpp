#!/bin/bash

set -v

tmpdir=tmpdir
mkdir -p $tmpdir

cp Tests/dcicpp.spc $tmpdir
cd $tmpdir
c=0
for problem in $(ls $MASTSIF/*.SIF)
do
  g=$(rundcicpp -D $problem -lgfortran -lgfortranbegin | grep Converged)
  if [ ! -z "$g" ]; then
    c=$(($c+1))
  fi
done

echo "Convergence count: $c"
[ $c -lt 700 ] && exit 1
