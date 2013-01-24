#!/bin/bash
list="PartialPenal project_dcp project_bfgs trustWorstdn trustConvexBox penal_trust penal_bfgs"
listsize=4
count=0
rm -rf directory*
while [ $count -le $[2**$listsize-1] ]
do
  #echo $[(count/4)%2]$[(count/2)%2]$[count%2]
  n1=$[count%2]
  n2=$[(count/2)%2]
  n3=$[(count/4)%2]
  n4=$[(count/8)%2]
  #n5=$[(count/16)%2]
  #n6=$[(count/32)%2]
  #n7=$[(count/64)%2]
  #echo $n7$n6$n5$n4$n3$n2$n1
  #number=$n7$n6$n5$n4$n3$n2$n1
  number=$n4$n3$n2$n1
  dirname=directory$number
  rm -f dcicpp.spc
  echo PartialPenal $n1 >> dcicpp.spc 
  echo project_dcp 0 >> dcicpp.spc
  #echo project_dcp $n2 >> dcicpp.spc 
  echo project_bfgs $n2 >> dcicpp.spc 
  echo trustWorstdn $n3 >> dcicpp.spc 
  echo trustConvexBox $n4 >> dcicpp.spc 
  echo penal_trust 0 >> dcicpp.spc
  echo penal_bfgs 1 >> dcicpp.spc
  echo ScaleVertical 1 >> dcicpp.spc
  #echo penal_trust $n6 >> dcicpp.spc 
  #echo penal_bfgs $n7 >> dcicpp.spc 
  mkdir -p $dirname
  cp runtest.sh $dirname
  cp dcicpp.spc $dirname
  (cd $dirname; ./runtest.sh ../TestLists/HS.list)
  #(cd $dirname; runcppcuter -p dcicpp -D HS35)
  #mv latex* $dirname
  count=$[count+1]
done
