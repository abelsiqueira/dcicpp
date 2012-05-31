#!/bin/bash
git log -1 --pretty=format:%H%n >> goodResults
wc -l latex* >> goodResults
echo ' ' >> goodResults
cat dcicpp.spc >> goodResults
