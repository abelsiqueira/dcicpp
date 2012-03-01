#/bin/sh
#for pname in $(cat littleTest)
for pname in $(cat TestLists/small.list)
do
  runcppcuter -p dcicpp -D $pname
done
