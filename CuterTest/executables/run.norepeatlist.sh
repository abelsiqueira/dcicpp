#/bin/sh
for pname in $(cat $1)
do
  grep_output="$(grep $pname latex_*)"
  [ "$grep_output" ] && echo "Skipping $pname" && continue
  echo "[$(date +%X)]Running $pname"
  runcppcuter -p dcicpp -D $pname
done
