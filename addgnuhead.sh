#!/bin/bash
for file in $(find . -name *cpp)
do
  cat GNU.head $file > $file.tmp
  mv $file.tmp $file
done
