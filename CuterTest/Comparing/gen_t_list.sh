#!/bin/bash

if [ $# -lt 1 ]
then
  echo "ERROR"
  exit
fi

out=$(echo $1 | sed 's/.ratio/.t_list/g')
echo $out
./gen_t_list.py $1 | sort -nu > $out
