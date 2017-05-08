# loop through all .h and .cc files
#for f in functions_unittest.cc; do
for f in `ls *.h *.cc`; do 
  # # double space on header guards
  # sed 's/#endif \/\//#endif  \/\//' $f > tttmp
  # mv tttmp $f

#  # remove extra whitespace at end of line
#  sed 's/ $//' $f > tttmp
#  mv tttmp $f

  # rename functions
  sed 's/myRanInit/ranInit/g' $f > tttmp
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
  mv tttmp $f

done
