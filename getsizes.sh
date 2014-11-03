#!/bin/bash

[ $# -lt 1 ] && echo "Need argument" && exit 1
file=$1
[ ! -f $file ] && echo "$file is not a file" && exit 1

[ -z "$(echo $file | grep "list$")" ] && echo "$file is not a .list file" && \
  exit 1

outfile=$(echo $file | sed 's/.list/.sizes/g')
rm -f $outfile

for p in $(cat $file)
do
  sifdecoder -D $p | awk -v p=$p -v v=0 -v c=0 \
    '/variables/ {v=v+$3};
     /constraints/ {c=c+$3};
     END{print p, v, c}' >> $outfile
done
