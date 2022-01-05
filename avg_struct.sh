traj=$1
structure=$2
base=$(basename $traj .xtc)
dir=$(dirname $traj)

gmx rmsf -f $file -s $structure -ox ${dir}/${base}_avg.pdb<<EOF
0
EOF

script_path=$(dirname "$0")

scwrl_path=$3

#run scwrl
python $script_path/run_scwrl.py -s ${dir}/${base}_avg.pdb --scwrl4 $scwrl_path

gmx pdb2gmx -f ${dir}/${base}_avg_scwrl.pdb -o ${dir}/${base}_avg_scwrl_prep.pdb -ignh<<EOF
6
1
EOF

python $script_path/finalEnergyMinim.py ${dir}/${base}_avg_scwrl_prep.pdb