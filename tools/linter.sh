# loop through all .h and .cc files
for f in `ls *.h *.cc`; do 
  # # double space on header guards
  # sed 's/#endif \/\//#endif  \/\//' $f > tttmp
  # mv tttmp $f

  # remove extra whitespace at end of line
  sed 's/ $//' $f > tttmp
  mv tttmp $f

done
