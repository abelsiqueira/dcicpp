#!/bin/bash

set -ev

tmpdir=tmpdir
mkdir -p $tmpdir

cp Tests/dcicpp.spc $tmpdir
cd $tmpdir
c=0
for n in $(seq 1 119)
do
  g=$(rundcicpp -D HS$n | grep Converged)
  if [ ! -z "$g" ]; then
    c=$[$c+1]
  fi
done

[ $c -lt 110 ] && exit 1
