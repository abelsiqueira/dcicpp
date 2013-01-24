#/bin/sh
for pname in $(ls latex*)
do
  grep unc $pname > Results/$pname'_unc'
  grep con $pname > Results/$pname'_con'
done
