#!/bin/bash

cat $1 | awk 'NR == 1 {max=$1} \
  { if ($1>max) max=$1} END {printf "%f", max}'
