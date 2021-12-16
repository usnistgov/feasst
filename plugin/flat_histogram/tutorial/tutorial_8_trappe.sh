#!/bin/bash

((num_hours=1))
#((num_hours=5*24))
num_procs=12
((num_procs_ext=2*$num_procs))

function launch_node {
# write output script to file
cat << _EOF_ > launch.cmd
#!/bin/bash
#SBATCH -n ${num_procs_ext}
#SBATCH -N 1
#SBATCH -t ${num_hours}:00:00
#SBATCH -o hostname_%j.out
#SBATCH -e hostname_%j.out
echo "Running on host \$(hostname)"
echo "Time is \$(date)"
echo "Directory is \$PWD"
echo "ID is \$SLURM_JOB_ID"

cd \$PWD

python tutorial_8_trappe.py \
  --task \$SLURM_ARRAY_TASK_ID \
  -p ~/feasst/forcefield/ethane.fstprt \
  -p ~/feasst/forcefield/ethene.fstprt \
  --num_hours $num_hours \
  --num_procs $num_procs

## for C++
#./tutorial \
#  --task \$SLURM_ARRAY_TASK_ID \
#  --particle0 ~/feasst/forcefield/ethane.fstprt \
#  --particle1 ~/feasst/forcefield/ethene.fstprt \
#  --num_hours $num_hours \
#  --num_procs $num_procs

## cylinder example
#./tutorial --task \$SLURM_ARRAY_TASK_ID \
#  --particle0 ~/feasst/forcefield/ethane.fstprt \
#  --particle1 ~/feasst/forcefield/ethene.fstprt \
#  --num_hours $num_hours \
#  --num_procs $num_procs \
#  --lx 50 \
#  --ly 50 \
#  --lz 50 \
#  --cyl_radius 25 \
#  --cyl_rcut 10 \
#  --max_particles 50

echo "Job is done"
echo "Time is \$(date)"
_EOF_

sbatch --array=0-15%1 launch.cmd
}

## for C++
#cmake .
#make

launch_node
