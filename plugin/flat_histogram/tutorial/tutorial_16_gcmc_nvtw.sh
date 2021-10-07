#!/bin/bash

((num_hours=5*24))
num_procs=32
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

python tutorial_16_gcmc_nvtw.py --temperature $temperature  --max_particles $max_particles --num_hours $num_hours --num_procs $num_procs --model $model --mu $mu

echo "Job is done"
echo "Time is \$(date)"
_EOF_

sbatch --array=0-3%1 launch.cmd
}

#model=sqw;temperature=0.7; max_particles=450; mu=-2.352321;
model=spce;temperature=300; max_particles=296; mu=-38;

launch_node
