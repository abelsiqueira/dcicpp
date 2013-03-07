#!/bin/bash
if [ $# -lt 1 ]
then
  echo "ERROR: Forgot the test name"
  exit 1
fi
function color_echo () {
  echo -e "\033[0;32m"$1"\033[0m"
}

color_echo "Calculating ratio"
./calc_performance_ratio.py $1 $2

color_echo "Generating list of values for t"
./gen_t_list.sh $1.ratio

color_echo "Removing symbolic value for maximum from t"
sed -i '1d' $1.t_list

color_echo "Obtaining maximum value for graph"
max=$(echo "scale=4; $(./find_max.py $1.t_list) + 1" | bc)

color_echo "Substituting symbolic value to maximum value"
sed -i "s/max/${max}/g" $1.ratio

color_echo "Generating m code for octave"
./gen_perf_func.py $1 > $1.m

color_echo "Executing octave to obtaing png"
echo `octave -q --eval "$1"`
