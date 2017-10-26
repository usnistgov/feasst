#!/usr/bin/env bash
SRC=$1
DES="$(basename $SRC .cc)".py
TMP=/tmp/tmp

grep -v "^#include \"" $SRC > $TMP
grep -v "^int main(" $TMP > $DES

sed --in-place 's/^  //' $DES
sed --in-place 's/;//g' $DES
sed --in-place 's/\/\//#/g' $DES
sed --in-place 's/feasst::\(.*\) \(.*\)(/\2 = feasst.\1(/g' $DES
sed --in-place 's/feasst::/feasst\./g' $DES
sed --in-place 's/vector<\(.*\)> \(.*\)(\(.*\))/\2 = list(range(\3))/g' $DES
sed --in-place 's/<</+/g' $DES
sed --in-place 's/const double //g' $DES
sed --in-place 's/const int //g' $DES
sed --in-place 's/\.str()//g' $DES
sed --in-place 's/\.c_str()//g' $DES
sed --in-place 's/exp(/math\.exp(/g' $DES
sed --in-place 's/&//g' $DES
sed --in-place 's/->/\./g' $DES
sed --in-place 's/}//g' $DES

cat <<- EOF > $TMP
import feasst
EOF

cat $DES >> $TMP
mv $TMP $DES

