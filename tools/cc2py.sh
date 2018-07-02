#!/usr/bin/env bash
# this script converts a cc file to python, typically used for tutorials.
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
sed --in-place 's/1e\(.[1-9]*\)/int(1e\1)/' $DES
sed --in-place 's/\/\*\*/"""/' $DES
sed --in-place 's/ \*\//"""\n\nimport feasst/' $DES
sed --in-place 's/double(/float(/' $DES
sed --in-place 's/pow(\(.*\), \(.*\))/(\1)**(\2)/' $DES
sed --in-place 's/^auto //' $DES
sed --in-place 's/{"\(.[a-z,A-Z]*\)", /"\1" : /' $DES
sed --in-place '/^std::stringstream/d' $DES
sed --in-place 's/ASSERT(\(.*\), "\(.*\)")/assert(\1) # \2/' $DES
sed --in-place 's/feasst\.str\(/str\(' $DES

#cat <<- EOF > $TMP
#import feasst
#EOF

#cat $DES >> $TMP
#mv $TMP $DES

