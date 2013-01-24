#!/bin/bash

# EDIT dcicpp.spc manually

MAXIT='200'
MAXTIME='60'

for file in $(ls TestLists/testes.[eiug]*)
do 
  ./clatex.sh
  outname=$(echo $file | sed 's/TestLists\///g')
  ./runlist.sh $file > $outname.out.$MAXIT.$MAXTIME
  wc -l latex* > $outname.wc.$MAXIT.$MAXTIME
  #for all results. 
  #for lfile in $(ls latex_*)
  #for fails only
  for lfile in $(ls latex_[^c]*)
  do
    echo $lfile
    ./take_first_arg.py $lfile
  done > $outname.fail.$MAXIT.$MAXTIME
done

