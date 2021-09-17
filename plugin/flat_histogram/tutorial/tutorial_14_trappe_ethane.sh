function launch_node {
# write output script to file
cat << _EOF_ > launch.cmd
#!/bin/bash
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -t 5-00:00
#SBATCH -o hostname_%j.out
#SBATCH -e hostname_%j.out
echo "Running on host \$(hostname)"
echo "Time is \$(date)"
echo "Directory is \$PWD"
echo "ID is \$SLURM_JOB_ID"

cd \$PWD

python tutorial_8_trappe.py -p ~/feasst/forcefield/ethane.fstprt --lx 30 --ly 30 --lz 30 --temperature 300 --num_procs 12 --max_particles 223 --beta_mu -7

echo "Job is done"
echo "Time is \$(date)"
_EOF_

sbatch launch.cmd
}

launch_node
