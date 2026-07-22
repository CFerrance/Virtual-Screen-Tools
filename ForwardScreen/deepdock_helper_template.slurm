#!/bin/bash
#SBATCH --job-name=deepdock_helper
#SBATCH --partition=CHOOSE_PARTITION
#SBATCH --nodelist=CHOOSE_NODE
#SBATCH --nodes=1
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=END,FAIL

echo Running on `hostname`
echo workdir $SLURM_SUBMIT_DIR
echo Start Time: $(date)

#Constant Paths
CLUSTER_YAML=slurm.yaml         #from GNINA directory
SPRINT=/SPRINT/CHECKPOINT/PATH  
GNINA=/GNINA/EXE/PATH


COMMAND=$1              #deepdock command
PDB_PROTEIN=$2          #receptor pdb file
PDB_LIGAND=$3           #ligand pdb file
SEQUENCE=$4             #receptor fasta file
CHEMICAL_LIBRARY=$5     #chemical library .parquet folder
OUTPUT_DIR=$6

#scratch drive folder to work in
SCRDIR=/net/node/scr/user

#if the scratch drive doesn't exist (it shouldn't) make it.
if [[ ! -e $SCRDIR ]]; then
        mkdir $SCRDIR
fi
chmod +rX $SCRDIR
echo scratch drive ${SCRDIR}

#Copy Inputs
cp $SLURM_SUBMIT_DIR/$PDB_PROTEIN ${SCRDIR}
cp $SLURM_SUBMIT_DIR/$PDB_LIGAND ${SCRDIR}
cp $SLURM_SUBMIT_DIR/$SEQUENCE ${SCRDIR}
cp $SLURM_SUBMIT_DIR/$CLUSTER_YAML ${SCRDIR}

cd ${SCRDIR}

#setup to copy files back to working dir on exit
trap "mv $OUTPUT_DIR $SLURM_SUBMIT_DIR" EXIT

python $SLURM_SUBMIT_DIR/deepdock.py $COMMAND \
    --input $CHEMICAL_LIBRARY \
    --dir ${SCRDIR}/$(basename $OUTPUT_DIR) \
    --cluster ${SCRDIR}/$(basename $CLUSTER_YAML) \
    --target_sequence ${SCRDIR}/$(basename $SEQUENCE) \
    --sprint_checkpoint $SPRINT \
    --receptor ${SCRDIR}/$(basename $PDB_PROTEIN) \
    --autobox_ligand ${SCRDIR}/$(basename $PDB_LIGAND) \
    --gnina_executable $GNINA \

echo Finish Time: $(date)