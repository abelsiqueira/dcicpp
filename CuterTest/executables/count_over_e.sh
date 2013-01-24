#!/bin/bash
total_val=0
total_total=0
for file in $(ls cond_*)
do
  val=$(grep "e+" $file | wc -l)
  total_val=$[$total_val + $val]
  total=$(cat $file | wc -l)
  total_total=$[$total_total + $total]
  if [ $total -gt 0 ]
  then
    porc=$(echo "scale=2; 100*$val/$total" | bc)
    echo -e "$file:   \t$val de $total   \t $porc%"
  fi
done
porc=$(echo "scale=2; 100*$total_val/$total_total" | bc)
echo -e "TOTAL:  \t\t$total_val de $total_total \t $porc%"
