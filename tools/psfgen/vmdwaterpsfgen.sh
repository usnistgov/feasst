l="25.89"  #gives 512 water molecules
l="12.19"  #gives 52 water molecules
#echo "solvate $l" > vmd
echo "
package require solvate
solvate -minmax {{0 0 0} {$l $l $l}}
#animate write psf wat.psf
quit
" > vmd

vmd -e vmd -dispdev none


