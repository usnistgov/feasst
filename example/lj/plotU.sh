filename2="$(basename $0 .sh)p${p}c${col}"

gnuplot << EOF
set terminal postscript eps enhanced color "Arial"
set encoding iso_8859_1
set size 0.65,0.55
set output "${filename2}.eps"
set xlabel "N"
set ylabel "U/N"
set key top right
p 'logp0' u 2:3,\
  'logp11' u 2:3,\
  'logp6' u 2:3,\
  'colMat' u 1:(\$3/\$1) w lines lc -1 lt -1
EOF

  #produce jpegs
  gs -sDEVICE=jpeg -dJPEGQ=100 -dNOPAUSE -dBATCH -dSAFER -r300 -sOutputFile=${filename2}.jpg ${filename2}.eps
  mogrify -trim -resize 1800x1600 ${filename2}.jpg
  display -resize 600x600 ${filename2}.jpg &


