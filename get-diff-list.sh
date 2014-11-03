#!/bin/bash

function usage {
  echo "$0 FILE1 FILE2"
}

function check_file {
  [ ! -f $1 ] && echo "$1 is not a file" && exit 1
}

[ $# -lt 2 ] && usage && exit 1
file1=$1
file2=$2

check_file $file1
check_file $file2

all=$(sort -u $file1 $file2)

rm -f both-fail.list only-$file1 only-$file2
for p in $all
do
  g1=$(grep $p $file1)
  g2=$(grep $p $file2)
  if [ ! -z "$g1" -a ! -z "$g2" ]; then
    echo $p >> both-fail.list
  elif [ ! -z "$g1" ]; then
    echo $p >> only-$file1
  elif [ ! -z "$g2" ]; then
    echo $p >> only-$file2
  fi
done

