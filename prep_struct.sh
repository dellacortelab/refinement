#script to prep structures for input into OpenMM using Gromacs. Give input structure path as first argument, out path as second argument. 
print_help () {
    echo -e "\nScript to use Gromacs to prep a structure for input into OpenMM. \n"
    echo "[-h]: print this help text and exit."
    echo "Required Arguments: "
    echo "[-s STRUCT_PATH]: path to structure to prepare."
    echo "[-o OUT_PATH]: path to write prepped structure out to."
    echo -e "\nSpecial Dependencies: Gromacs\n"
    exit 1
}

if [ $# -eq 0 ]; then
    print_help
fi

has_struct=false
has_out=false

IFS=' ' read -r -a allArgs <<< "$@"
for ((idx=0;idx<$#; idx++)); do
    if [ ${allArgs[${idx}]} == "-h" ]; then
        print_help
    elif [ ${allArgs[${idx}]} == "-s" ]; then
        structure=${allArgs[$(($idx+1))]}
        has_struct=true
    elif [ ${allArgs[${idx}]} == "-o" ]; then
        out_path=${allArgs[$(($idx+1))]}
        has_out=true
    fi
done

if [ "$has_struct" = false ]; then
    echo "Missing required argument -s"
    exit
elif [ "$has_out" = false ]; then
    echo "Missing required argument -o"
    exit
fi

gmx pdb2gmx -f $structure -o $out_path -ignh<<EOF
6
1
EOF