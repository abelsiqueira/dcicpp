#/bin/sh
for pname in $(cat $1)
do
  runcppcuter -p dcicpp -D $pname
done
