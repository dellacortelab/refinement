#script to prep structures for input into OpenMM using Gromacs. Give input structure path as first argument, out path as second argument. 

gmx pdb2gmx -f $1 -o $2 -ignh<<EOF
6
1
EOF