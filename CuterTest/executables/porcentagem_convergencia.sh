#!/bin/bash
target_dir=$1
( wc -l $target_dir/latex_* latex_* | grep conv  | awk '{s = s + $1}END{print s}'; 
  wc -l $target_dir/latex_* latex_* | grep total | awk '{s = s + $1}END{print s}' ) 
  | while read i; do echo -n "$i "; done | awk '{print $1/$2}'
