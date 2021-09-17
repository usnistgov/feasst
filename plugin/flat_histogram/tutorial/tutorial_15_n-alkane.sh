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

python tutorial_15_n-alkane.py --particle $data --temperature $temperature  --max_particles $max_particles --lx $box_length --ly $box_length --lz $box_length --beta_mu $beta_mu --cutoff $cutoff --task \$SLURM_ARRAY_TASK_ID --num_hours $num_hours --num_procs $num_procs

echo "Job is done"
echo "Time is \$(date)"
_EOF_

sbatch --array=0-3%1 launch.cmd
}

data=ethane.fstprt; temperature=300; max_particles=225; box_length=30; beta_mu=-7.1126; cutoff=15 ## https://mmlapps.nist.gov/srs/ETHANE/ethane_sat.htm
#data=~/feasst/forcefield/propane.fstprt; temperature=312; max_particles=200; box_length=30; beta_mu=-7.1126; cutoff=14;
#data=~/feasst/forcefield/propane.fstprt; temperature=344; max_particles=180; box_length=28; beta_mu=-7; cutoff=14;
#data=~/feasst/forcefield/n-butane.fstprt; temperature=350; max_particles=575; box_length=45; beta_mu=-6; cutoff=12; # https://www.nist.gov/mml/csd/chemical-informatics-research-group/sat-tmmc-liquid-vapor-coexistence-properties-trappe-ua-n
#data=~/feasst/forcefield/n-octane.fstprt; temperature=560; max_particles=180; box_length=35; beta_mu=-6; cutoff=15; # https://mmlapps.nist.gov/srs/OCTANE/octane_sat.htm
#data=~/feasst/forcefield/n-decane.fstprt; temperature=560; max_particles=120; beta_mu=-7; box_length=36;

launch_node
