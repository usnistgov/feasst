for f in *; do
  #sed 's/get_vector(/coord(/g' $f > ttmp; mv ttmp $f
  #sed 's/get_//g' $f > ttmp; mv ttmp $f
  #sed 's/ "include/ "core\/include/' $f > ttmp; mv ttmp $f
  #sed 's/ "src/ "core\/src/' $f > ttmp; mv ttmp $f
  #sed 's/ "core\/scr/ "core\/src/' $f > ttmp; mv ttmp $f
  #sed 's/ EWALD_INCLUDE_/ FEASST_EWALD_/' $f > ttmp; mv ttmp $f
  #sed 's///g' $f > ttmp; mv ttmp $f
  sed 's/running_energy/current_energy/g' $f > ttmp; mv ttmp $f
  #sed 's/DomainCuboid/Domain/g' $f > ttmp; mv ttmp $f
done
#rename 's/\.cc/\.cpp/' $(find . -type f)
#rename -n 's/\.cc/\.cpp/' $(find . -type f)

#for f in *.cc; do
#  git mv "$f" "$(basename "$f" .cc).cpp"
#done
