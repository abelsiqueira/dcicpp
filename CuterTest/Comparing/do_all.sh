#!/bin/bash
# Primeiro eh preciso pegar os resultados e criar o ratio
./calc_performance_ratio.py $1
# Agora preciso calcular t
./gen_t_list.sh $1.ratio
sed -i '1d' $1.t_list
max=$(echo "scale=4; $(./find_max.py $1.t_list) + 1" | bc)
sed -i "s/max/${max}/g" $1.ratio
# Agora gera o .m
./gen_perf_func.py $1 > $1.m
# Executa o .m 
echo `octave -q --eval "$1"`
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$1.pdf $1.ps
