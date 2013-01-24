#/bin/sh
rm -f latex*
for pname in $(cat TestLists/small.list)
do
  runcppcuter -p dcicpp -D $pname
done
