#!/bin/bash

# EDIT dcicpp.spc manually

target_dir=all.small.2013.02.19.saving

for line in $(grep Skipping $target_dir/*.out | sed 's/ /:/g')
do
  outname=$(echo $line | awk -F ':' '{print $1}' | sed 's/.out//g')
  problem_name=$(echo $line | awk -F ':' '{print $3}')
  runcppcuter -p dcicpp -D $problem_name >> $outname.out
  wc -l latex* >> $outname.extra.wc
  for lfile in $(ls latex_*)
  do
    cat $lfile >> $target_dir/$lfile.extra
  done
  ./executables/clatex.sh
done


exit
for list in $test_list
do 
  outname=$(echo $list | sed "s:TestLists/small/classification.:$target_dir/all.small.:g")
  [ -f $outname.wc ] && echo "skipping $outname" && continue
  echo "Running $outname"
  ./executables/run.norepeatlist.sh $list >> $outname.out
  wc -l latex* > $outname.wc
  for lfile in $(ls latex_*)
  do
    cat $lfile >> $target_dir/$lfile
    echo $lfile
    ./executables/take_first_arg.py $lfile
  done > $outname.list.of.problems
  echo "error solving:" >> $outname.list.of.problems
  for problem in $(cat $list)
  do
    grep_output="$(grep $problem latex_*)"
    if [ -z "$grep_output" ]
    then
      echo $problem
    fi
  done >> $outname.list.of.problems
  ./executables/clatex.sh
done

