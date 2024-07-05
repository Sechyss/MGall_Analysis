#!/bin/bash -l

#SBATCH --partition=defq
#SBATCH --job-name=BBmap
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks=16                         # Run a single task (increase this value for parallelisation across CPUs)
#SBATCH --account=c.bonneaud                # The accounting code - usually named after the PI for a project
#SBATCH --constraint=IB                     # Specify features required by the job
#SBATCH --mem=23000                         # Memory per node specification is in MB. It is optional.
#SBATCH --output=bbmap_Lucy-%j.out          # Standard output and error log
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=at991@exeter.ac.uk      # E-Mail address of the user that needs to be notified.
##SBATCH --requeue                          # Specifies that the job will be requeued after a node failure.
                                            # The default is that the job will not be requeued.

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

module load BBTools/38.39
module load samtools/1.10

cd /nobackup/beegfs/workspace/at991/Data/Lucy_reads/ || exit
bbmap.sh ref=/nobackup/beegfs/workspace/at991/Data/R.fna

# List all files in the directory and process each one
for file in /nobackup/beegfs/workspace/at991/Data/Lucy_reads/*_R1.fastq.gz; do
  # Extract the filename from the full path
  filename=$(basename "$file")

  # Cut the filename after the string 'sickle'
  newname="${filename%%_sickle*}"

  # Define the output directory
  output_dir="/nobackup/beegfs/workspace/at991/Data/Mapped_output/${newname}.sam"

  bbwrap.sh build=1 in1=/nobackup/beegfs/workspace/at991/Data/Lucy_reads/"$newname"sickle_R1.fastq.gz,/nobackup/beegfs/workspace/at991/Data/Lucy_reads/"$newname"sickle_single.fastq.gz in2=/nobackup/beegfs/workspace/at991/Data/Lucy_reads/"$newname"sickle_R2.fastq.gz,null -out="$output_dir"
done

cd /nobackup/beegfs/workspace/at991/Data/Mapped_output/ || exit

for file in /nobackup/beegfs/workspace/at991/Data/Mapped_output/*sam; do
  # Extract the filename from the full path
  filename=$(basename "$file")

  # Cut the filename after the string 'nophi'
  newname="${filename%%_nophi_*}"

  # Define the output directory
  output_dir="/nobackup/beegfs/workspace/at991/Data/Mapped_output/${newname}.bam"

  samtools sort -l 8 -o "$output_dir" -O BAM --threads 16 /nobackup/beegfs/workspace/at991/Data/Mapped_output/"$newname"_nophi_.sam

done