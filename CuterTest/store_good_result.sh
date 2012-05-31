#!/bin/bash
git show --pretty=format:%H%n >> goodResults
wc -l latex* >> goodResults
echo ' ' >> goodResults
cat dcicpp.spc >> goodResults
