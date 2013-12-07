#!/bin/bash

[ $# -lt 1 ] && echo "Need list" && exit 1
target_dir=cuter.$(date +"%Y.%m.%d_%H.%M")
list=$1

mkdir -p $target_dir
cp $list $target_dir/
cp dcicpp.spc $target_dir/

for problem in $(cat $list)
do
  runcppcuter -p dcicpp -D $problem > $target_dir/$problem.out
done
