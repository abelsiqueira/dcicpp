#!/bin/bash
outname='cond.inf'
#rm -f $outname
ls matrix_jacob*.m | while read file
do
  file=$(echo $file | sed 's/\.m//g')
  cond=$(octave -q --eval "$file; cond(A)")
  cond=$(echo $cond | sed 's/ans = //g')
  name=$(echo $file | sed 's/matrix_jacob//g')
  echo "Problem $name - k(A) = $cond" >> $outname
done
