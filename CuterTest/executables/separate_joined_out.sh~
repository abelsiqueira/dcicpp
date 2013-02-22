#!/bin/bash

[ $# -lt 1 ] && echo "Need joined.out clean" && exit 1

prefix=$(echo $1 | sed 's/.out//g')

name_file=$prefix.name
exec_time_file=$prefix.exec.time
setup_time_file=$prefix.setup.time
exit_file=$prefix.exit

grep sifdecode $1 | awk '{print $2}' > $name_file
grep "Execution time" $1 | awk '{print $3}' > $exec_time_file
grep "Setup time" $1 | awk '{print $3}' > $setup_time_file
grep "EXIT" $1 | awk '{print $5}' > $exit_file

wc -l $name_file
wc -l $exec_time_file
wc -l $setup_time_file
wc -l $exit_file
