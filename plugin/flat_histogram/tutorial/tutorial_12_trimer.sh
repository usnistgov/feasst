#!/bin/bash

# Launch a FEASST simulation on an HPC cluster that uses a SLURM queue.
# The SLURM array automatically restarts the simulation by the "task".
# The first job in the array has a task of 0, which signals the FEASST
# simulation to initialize the simulation.
# The remaining jobs in the array supply task > 0, which signals simulation
# restart using checkpoint files instead.

((num_hours=8))
num_procs=12
((num_procs_ext=2*$num_procs))

# Write a SLURM script to file and queue it.
function launch_node {
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
echo "TASK is \$SLURM_ARRAY_TASK_ID"

cd \$PWD
python tutorial_12_trimer.py --task \$SLURM_ARRAY_TASK_ID --num_procs ${num_procs} --num_hours $num_hours

if [ \$? == 0 ]; then
  echo "Job is done"
  scancel \$SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi

echo "Time is \$(date)"
_EOF_

sbatch --array=0-10%1 launch.cmd
}

launch_node
