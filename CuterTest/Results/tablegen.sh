#/bin/sh
rm -f result.table
touch result.table

wc -l latex_* > result.table
