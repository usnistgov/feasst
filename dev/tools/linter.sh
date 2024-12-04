#!/bin/bash

#for f in *; do
for f in $(find . -name '*.cpp' -o -name '*.h' -o -name '*.py' -o -name '*.dot' -o -name '*.ipynb'); do
  #sed 's/get_vector(/coord(/g' $f > ttmp; mv ttmp $f
  #sed 's/get_//g' $f > ttmp; mv ttmp $f
  #sed 's/ "include/ "core\/include/' $f > ttmp; mv ttmp $f
  #sed 's/ "src/ "core\/src/' $f > ttmp; mv ttmp $f
  #sed 's/ "core\/scr/ "core\/src/' $f > ttmp; mv ttmp $f
  #sed 's/ EWALD_INCLUDE_/ FEASST_EWALD_/' $f > ttmp; mv ttmp $f
  #sed 's///g' $f > ttmp; mv ttmp $f
  #sed 's/_CORE_/_MONTE_CARLO_/g' $f > ttmp; mv ttmp $f
  #sed 's/_CORE_/_SYSTEM_/g' $f > ttmp; mv ttmp $f
  #sed 's/core\/include/system\/include/g' $f > ttmp; mv ttmp $f
  #sed 's/core\/test/system\/test/g' $f > ttmp; mv ttmp $f
#  sed 's/_STEPPERS_/_MONTE_CARLO_/g' $f > ttmp; mv ttmp $f
  #sed 's/utils_file\.h/file.h/g' $f > ttmp; mv ttmp $f
  #sed 's/\&\&/and/g' $f > ttmp; mv ttmp $f
  #rm $f/ttmp
  #sed 's/(Potential(/(MakePotential(/g' $f > ttmp; mv ttmp $f
  #sed 's/steps_per/trials_per/g' $f > ttmp; mv ttmp $f
  sed 's/ num_iterations_to_complete / cycles_to_complete /g' $f > ttmp; mv ttmp $f
  #sed 's/{"num_trials_per_iteration"/{"trials_per_cycle"/g' $f > ttmp; mv ttmp $f
  #sed 's/{{"until_criteria_complete", "true"}}/{{"until", "complete"}}/g' $f > ttmp; mv ttmp $f
  #sed 's/Run until_criteria_complete true/Run until complete/g' $f > ttmp; mv ttmp $f
  #sed 's/AddReference/ConvertToRefPotential/g' $f > ttmp; mv ttmp $f
  #sed 's/steps_since/trials_since/g' $f > ttmp; mv ttmp $f
  #sed 's/epsilon/Epsilon/g' $f > ttmp; mv ttmp $f
  #sed 's/cutoff/CutOff/g' $f > ttmp; mv ttmp $f
  #sed 's/charge/Charge/g' $f > ttmp; mv ttmp $f
  #sed 's/patch_angle/PatchAngle/g' $f > ttmp; mv ttmp $f
  #sed 's/director/Director/g' $f > ttmp; mv ttmp $f
  #sed 's/criteria_metropolis\.h/metropolis\.h/g' $f > ttmp; mv ttmp $f
  #sed 's/ $//g' $f > ttmp; mv ttmp $f
done
#rename 's/\.cc/\.cpp/' $(find . -type f)
#rename -n 's/\.cc/\.cpp/' $(find . -type f)

#for f in *.cc; do
#  git mv "$f" "$(basename "$f" .cc).cpp"
#done
