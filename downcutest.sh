#!/bin/bash
cmd="svn checkout -q --username anonymous"
url="http://ccpforge.cse.rl.ac.uk/svn/cutest/"
for name in archdefs sifdecode cutest sif
do
  echo "Downloading $name"
  $cmd $url/$name/trunk ./$name
done
