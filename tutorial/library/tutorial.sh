#!/bin/bash

# Launch a FEASST simulation on an HPC cluster that uses a SLURM queue.

# Write a SLURM script to file and queue it.
function launch_node {
cat << _EOF_ > launch.cmd
#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -o hostname_%j.out
#SBATCH -e hostname_%j.out
echo "Running on host \$(hostname)"
echo "Time is \$(date)"
echo "Directory is \$PWD"
echo "ID is \$SLURM_JOB_ID"

cd \$PWD
python tutorial.py --trials $trials
#./tutorial --trials $trials

if [ \$? == 0 ]; then
  echo "Job is done"
  scancel \$SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi

echo "Time is \$(date)"
_EOF_

sbatch launch.cmd
}

## If C++, compile the tutorial.
#mkdir build
#cd build
#cmake ..
#make

trials=15000000

launch_node
