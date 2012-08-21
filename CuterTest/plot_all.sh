#!/bin/bash
mkdir -p Plots
ls plot_matlab*.m | while read file
do
  file=$(echo $file | sed 's/\.m//g')
  octave -q --eval \
    "$file; plot(ngp); hold on; plot(normc/max(normc),'r'); print('$file.ps','-color')"
  gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$file.pdf $file.ps
  mv $file.pdf Plots/
  rm -f $file.ps
done

