#!/bin/bash

# Compare the results with TraPPE CO2: http://trappe.oit.umn.edu/
# Comment on any observed differences. Possible causes include: Gibbs vs FH GCMC, system size, Ewald, etc

# Write a SLURM script to file and queue it.
function launch_node {
cat << _EOF_ > launch.cmd
#!/bin/bash
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -o hostname_%j.out
#SBATCH -e hostname_%j.out
echo "Running on host \$(hostname)"
echo "Time is \$(date)"
echo "Directory is \$PWD"
echo "ID is \$SLURM_JOB_ID"
echo "TASK is \$SLURM_ARRAY_TASK_ID"

cd \$PWD
python tutorial_7_spce_parallel.py --task \$SLURM_ARRAY_TASK_ID --num_procs 12 --num_hours 24 --molecule ~/feasst/forcefield/co2.fstprt --temperature 290 --max_molecules 300 --cubic_box_length 28

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
