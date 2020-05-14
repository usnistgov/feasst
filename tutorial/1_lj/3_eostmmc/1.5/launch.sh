# Custom script to launch jobs
output="$(basename $0 .sh)"
t=1.5; x=725; boxl=10; lnz="-1.568214"
for p in 0 1 2 3; do
  dir=t${t}x${x}l${boxl}z${lnz}p$p
  # write slurm script to file
cat << _EOF_ > ${output}.cmd
#!/bin/bash
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -o hostname_%j.out
#SBATCH -e hostname_%j.out
mkdir -p $dir/
cd $dir
mkdir -p tmp
echo "Time is \$(date)"
export OMP_NUM_THREADS=12
../../../../../tools/run.sh ../../test.py -t $t -z $lnz -x $x --openMP
echo "Time is \$(date)"
_EOF_

  sbatch ${output}.cmd
done
