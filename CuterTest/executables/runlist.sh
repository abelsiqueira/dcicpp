#!/bin/bash

[ $# -lt 1 ] && echo "Need list" && exit 1
target_dir=cutest.$(date +"%Y.%m.%d_%H.%M")
list=$1

mkdir -p $target_dir
cp $list $target_dir/
cp dcicpp.spc $target_dir/

cat >> $target_dir/INFORMATION << EOF
git localization
hash: $(git gethash)
status:
$(git status)
diff
$(git diff)
EOF

list=$(cat $list)
cd $target_dir
for problem in $list
do
  rundcicpp -p dcicpp -D $problem > $problem.out
done

echo -n "Problem & nvar & ncon " > table_latex
echo -n "& f & |gp| & |c|" >> table_latex
echo -n "& iter & time & type & bfgs" >> table_latex
echo "\\\\ \\hline" >> table_latex
sort latex_*  >> table_latex
