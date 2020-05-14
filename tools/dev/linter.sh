#!/usr/bin/env bash
# loop through all .h and .cc files
#for f in functions_unittest.cc; do
#for f in `ls *.i`; do

#find ../tutorial/ -name "test*" | xargs sed -i "s/’/'/"
#find . -type f | xargs sed -i "s/ $//g"

for f in `ls *`; do
#for f in `ls *.h *.cc *.i`; do
  # # double space on header guards
  # sed 's/#endif \/\//#endif  \/\//' $f > tttmp
  # mv tttmp $f

  # remove extra whitespace at end of line
  #sed --in-place 's/ $//g' $f

  # remove lines
  #sed --in-place 's/pair_lj_multi\.h/pair_lj\.h/g' $f
  #sed --in-place 's/PairLJMulti/PairLJ/g' $f

  # rename functions
  #sed --in-place '/#ifdef FEASST_NAMESPACE_/d' $f
  #sed --in-place '/#endif  \/\/ FEASST_NAMESPACE_/d' $f
  sed --in-place 's/field\.h/pair_field\.h/g' $f
  #sed --in-place 's/field_/pair_field_/g' $f
  #sed --in-place 's/FIELD_/PAIR_FIELD_/g' $f
  #sed --in-place 's/Field/PairField/g' $f
  # sed --in-place 's/print(/write(/g' $f
  #sed -i "s/’/'/g" $f
  #sed --in-place 's/void initEnergy;/void initEnergy();/g' $f

  #sed --in-place 's/namespace feasst {/\#ifdef FEASST_NAMESPACE_\nnamespace feasst {\n\#endif  \/\/ FEASST_NAMESPACE_/g' $f
  #sed --in-place 's/\}  \/\/ namespace feasst/\#ifdef FEASST_NAMESPACE_\n\}  \/\/ namespace feasst\n\#endif  \/\/ FEASST_NAMESPACE_/g' $f

  #sed 's/round(/feasstRound(/g' $f > tttmp
  #sed 's/myRanInit/ranInit/g' $f > tttmp
  #sed 's/myProd/feasst::product/g' $f > tttmp
  #sed 's/myMaxElement/feasst::maxElement/g' $f > tttmp
  #sed 's/myRound/feasst::round/g' $f > tttmp
  #sed 's/myAcos/feasst::arccos/g' $f > tttmp
  #sed 's/myTrim/feasst::trim/g' $f > tttmp
  #sed 's/myFileExists/fileExists/g' $f > tttmp
  #sed 's/myVecAv/vecAverage/g' $f > tttmp
  #sed 's/myMatMul/matMul/g' $f > tttmp
  #sed 's/myVecDotProd/vecDotProd/g' $f > tttmp
  #sed 's/myMatVecMul/matVecMul/' $f > tttmp
  #sed 's/myFill/feasst::fill/' $f > tttmp
  #mv tttmp $f

done
