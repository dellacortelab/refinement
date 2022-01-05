#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=10G   # memory per CPU core
#SBATCH -J "refinement"   # job name
#SBATCH --mail-user=con.morris09@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -C 'pascal'
# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

if [ $# -eq 0 ]; then
    echo "Do not submit this script on its own. Use runRefinement.sh"
    exit
fi


fullPath=$(pwd)
username=$(echo $fullPath | cut -d/ -f3)

refine=$2
restrain=$3
length=$4
gasPhase=$5
partialRestraint=$6
partialIdx=$7
hasGap=$8
resids=$9

module purge
module load cuda/10.1
if [ "$refine" = true ]; then
    module load gromacs/2019.4
fi

inputPDB="$1"
pdbBasename=$(basename $inputPDB .pdb)
#Allows for input w/ or w/out .pdb at end of filename
inputPDB="${pdbBasename}.pdb"
if [ "${refine}" = false ]; then
    targetname=$(echo "$inputPDB" | cut -d_ -f1)
    inPath="/fslhome/${username}/fsl_groups/fslg_dellacortelab/compute/CASP14/PRODUCTIVE_2020/reconstruction/${targetname}/postrelax/${inputPDB}"
    outPath="/fslhome/${username}/fsl_groups/fslg_dellacortelab/compute/CASP14/PRODUCTIVE_2020/refinement/${pdbBasename}/"
else
    targetname=$(basename ${pdbBasename} _prep)
    inPath="/fslhome/${username}/fsl_groups/fslg_dellacortelab/compute/CASP14/PRODUCTIVE_2020/refinement/${targetname}/${inputPDB}"
    outPath="/fslhome/${username}/fsl_groups/fslg_dellacortelab/compute/CASP14/PRODUCTIVE_2020/refinement/${targetname}/"
fi


echo $inputPDB

if [ "${refine}" = false ]; then
    mkdir $outPath
fi

if [ ${SLURM_ARRAY_TASK_ID} -eq 0 ] && [ "${refine}" = false ]; then
    cp $inPath $outPath
fi

scriptPath="/fslhome/${username}/fsl_groups/fslg_dellacortelab/compute/CASP14/refinement/scripts/parallel_refinement.py"

#requires OpenMM
echo "Starting Refinement Script"
cmd="python -u ${scriptPath} -i ${inPath} -o ${outPath} --iter ${SLURM_ARRAY_TASK_ID} --length $length"
if [ "$restrain" = false ]; then
    cmd="$cmd --norestraint"
    echo "Removing restraints"
fi

if [ "$gasPhase" = true ]; then
    cmd="$cmd --gas"
fi

if [ "$partialRestraint" = true ]; then
    cmd="$cmd -pr $partialIdx"
fi

if [ "$hasGap" = true ]; then
    cmd="$cmd -dr $resids"
fi

$cmd
    

