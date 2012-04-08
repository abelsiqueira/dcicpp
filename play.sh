#!/bin/bash
for i in $(cat dcicpp.spc)
do
  echo "case en_$i: aux >> $i; break;"
done
