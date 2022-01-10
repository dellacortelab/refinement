print_help () {
    echo -e "\n This script creates an averaged structure of a trajectory, then runs SCWRL4 to optimize sidechains and a restrained energy minimization on the averaged structure.\n"
    echo "[-h]: print this help text and exit."
    echo "Required Arguments: "
    echo "[-t TRAJ_PATH]: path to trajectory to create an averaged structure out of."
    echo "[-s STRUCT_PATH]: path to structure PDB file for trajectory."
    echo "[--scwrl4 SCWRL_PATH]: path to SCWRL4 executable."
    echo -e "\nSpecial Dependencies: Gromacs, SCWRL4, OpenMM\n"
    exit 1
}

if [ $# -eq 0 ]; then
    print_help
fi

has_traj=false
has_struct=false
has_scwrl=false

IFS=' ' read -r -a allArgs <<< "$@"
for ((idx=0;idx<$#; idx++)); do
    if [ ${allArgs[${idx}]} == "-h" ]; then
        print_help
    elif [ ${allArgs[${idx}]} == "-t" ]; then
        traj=${allArgs[$(($idx+1))]}
        has_traj=true
    elif [ ${allArgs[${idx}]} == "-s" ]; then
        structure=${allArgs[$(($idx+1))]}
        has_struct=true
    elif [ ${allArgs[${idx}]} == "--scwrl4" ]; then
        scwrl_path=${allArgs[$(($idx+1))]}
        has_scwrl=true
    fi
done

if [ "$has_traj" = false ]; then
    echo "Missing required argument -t"
    exit
elif [ "$has_struct" = false ]; then
    echo "Missing required argument -s"
    exit
elif [ "$has_scwrl" = false ]; then
    echo "Missing required argument --scwrl4"
    exit
fi

base=$(basename $traj .xtc)
dir=$(dirname $traj)

gmx rmsf -f $traj -s $structure -ox ${dir}/${base}_avg.pdb<<EOF
0
EOF

if [ $? -ne 0 ]; then
    echo 'gmx rmsf command failed to generate an averaged structure'
    exit
fi

script_path=$(dirname "$0")

#run scwrl
python $script_path/run_scwrl.py -s ${dir}/${base}_avg.pdb --scwrl4 $scwrl_path

if [ $? -ne 0 ]; then
    echo 'Failed to run scwrl successfully'
    exit
fi

gmx pdb2gmx -f ${dir}/${base}_avg_scwrl.pdb -o ${dir}/${base}_avg_scwrl_prep.pdb -ignh<<EOF
6
1
EOF

python $script_path/finalEnergyMinim.py ${dir}/${base}_avg_scwrl_prep.pdb

if [ $? -ne 0 ]; then
    echo 'Failed to energy minimize final structure'
    exit
fi

cp ${dir}/${base}_avg_scwrl_prep_minim.pdb ${dir}/${base}_final.pdb
if [ $? -eq 0 ]; then
    echo "Final refined structure saved to ${dir}/${base}_final.pdb"
fi