#!/bin/bash -l

#SBATCH --partition=defq
#SBATCH --job-name=iqtree
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks=20                         # Run a single task (increase this value for parallelisation across CPUs)
#SBATCH --account=c.bonneaud                # The accounting code - usually named after the PI for a project
#SBATCH --constraint=IB                     # Specify features required by the job
#SBATCH --mem=23000                         # Memory per node specification is in MB. It is optional.
#SBATCH --output=iqtree_result-%j.out       # Standard output and error log
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=at991@exeter.ac.uk      # E-Mail address of the user that needs to be notified.
##SBATCH --requeue                          # Specifies that the job will be requeued after a node failure.
                                            # The default is that the job will not be requeued.

module load ks575/IQ-tree/v2.0.3

cd /nobackup/beegfs/workspace/at991/Data/tree_alignments/ || exit

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR


for file in /nobackup/beegfs/workspace/at991/Data/tree_alignments/*.fasta; do
  # Extract the filename from the full path
  filename=$(basename "$file")
  newname="${filename%%.fasta*}"
  iqtree -s file --prefix "$newname" -m GTR+F+I+G4 -bb 1000 -t 20
done

