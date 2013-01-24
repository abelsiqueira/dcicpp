#!/bin/bash
listsize=8
count=0
while [ $count -le $[2**$listsize-1] ]
do
  #echo $[(count/4)%2]$[(count/2)%2]$[count%2]
  n1=$[count%2]
  n2=$[(count/2)%2]
  n3=$[(count/4)%2]
  n4=$[(count/8)%2]
  n5=$[(count/16)%2]
  n6=$[(count/32)%2]
  n7=$[(count/64)%2]
  n8=$[(count/128)%2]
  #echo $n7$n6$n5$n4$n3$n2$n1
  number=$n8$n7$n6$n5$n4$n3$n2$n1
  #number=$n3$n2$n1
  #dirname=allpossibilities$number
  rm -f dcicpp.spc
  echo PartialPenal $n1 >> dcicpp.spc 
  echo project_dcp $n2 >> dcicpp.spc 
  echo project_bfgs $n3 >> dcicpp.spc 
  echo trustWorstdn $n4 >> dcicpp.spc 
  echo trustConvexBox $n5 >> dcicpp.spc 
  echo penal_trust $n6 >> dcicpp.spc 
  echo penal_bfgs $n7 >> dcicpp.spc 
  echo ScaleVertical $n8 >> dcicpp.spc
  #mkdir -p $dirname
  #cp runtest.sh $dirname
  #cp dcicpp.spc $dirname
  #(cd $dirname; ./runtest.sh ../TestLists/HS.list)
  #(cd $dirname; runcppcuter -p dcicpp -D HS35)
  #mv latex* $dirname
  runcppcuter -p dcicpp -D $1 > $1.$number
  count=$[count+1]
done
