#!/bin/bash
#SBATCH --job-name=deepdock_helper
#SBATCH --partition=dept_cpu,dept_gpu
#SBATCH --nodelist=NODE
#SBATCH --nodes=1

echo Running on `hostname`
echo workdir $SLURM_SUBMIT_DIR
echo Start Time: $(date)

eval "$(conda shell.bash hook)"
conda activate deepdock_env

## CONSTANTS
## molPort library: /net/galaxy/home/koes/dkoes/data/molport.smiles
## enamine library: /net/g023/scr/dkoes/REAL/REAL.parquet
## cluster yaml: local slurm.yaml
## gnina executable: /net/galaxy/home/koes/dkoes/local/bin/gnina
## Sprint checkpoint: /net/galaxy/home/koes/dkoes/git/panspecies-dti/checkpoints/sprint.ckpt

COMMAND=$1              #deepdock command
MAX_WORKERS=$2          #max workers

DEEPDOCK=$HOME/path/to/deepdock.py
RECEPTOR=$HOME/path/to/receptor
LIGAND=$HOME/path/to/ligand
SEQUENCE=$HOME/path/to/sequence
LIBRARY=path/to/library             #replace with library from constants
CLUSTER_YAML=$HOME/path/to/slurm.yaml
OUT=deepdock_results

SCRDIR=/net/NODE/scr/USER           #scratch drive path, replace with node + folder name


#if the scratch drive doesn't exist (it shouldn't) make it.
if [[ ! -e $SCRDIR ]]; then
        mkdir $SCRDIR
fi
chmod +rX $SCRDIR
echo scratch drive ${SCRDIR}

#Copy Inputs to scratch directory for faster access
cp $RECEPTOR ${SCRDIR}
cp $LIGAND ${SCRDIR}
cp $SEQUENCE ${SCRDIR}
cp $CLUSTER_YAML ${SCRDIR}

cd ${SCRDIR}

#setup to copy files back to working dir on exit
trap "mv $OUTPUT_DIR $SLURM_SUBMIT_DIR" EXIT

python $DEEPDOCK $COMMAND \
    --input $LIBRARY \
    --dir ${SCRDIR}/$(basename $OUT) \
    --cluster ${SCRDIR}/$(basename $CLUSTER_YAML) \
    --target_sequence ${SCRDIR}/$(basename $SEQUENCE) \
    --sprint_checkpoint /net/galaxy/home/koes/dkoes/git/panspecies-dti/checkpoints/sprint.ckpt \
    --receptor ${SCRDIR}/$(basename $PDB_PROTEIN) \
    --autobox_ligand ${SCRDIR}/$(basename $PDB_LIGAND) \
    --gnina_executable /net/galaxy/home/koes/dkoes/local/bin/gnina \
    --max_workers=$MAX_WORKERS \
    --iolimit 100 \
    --smiles-only

echo Finish Time: $(date)

