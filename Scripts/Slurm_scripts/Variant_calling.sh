#!/bin/bash -l

#SBATCH --partition=defq
#SBATCH --job-name=VCF
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks=16                         # Run a single task (increase this value for parallelisation across CPUs)
#SBATCH --account=c.bonneaud                # The accounting code - usually named after the PI for a project
#SBATCH --constraint=IB                     # Specify features required by the job
#SBATCH --mem=23000                         # Memory per node specification is in MB. It is optional.
#SBATCH --output=vcf_Lucy-%j.out            # Standard output and error log
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=at991@exeter.ac.uk      # E-Mail address of the user that needs to be notified.
##SBATCH --requeue                          # Specifies that the job will be requeued after a node failure.
                                            # The default is that the job will not be requeued.

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

module load bcftools/1.9

cd /nobackup/beegfs/workspace/at991/Data/Mapped_output_WI01_2001_043_13_2P/ || exit

for file in /nobackup/beegfs/workspace/at991/Data/Mapped_output_WI01_2001_043_13_2P/*bam; do
  # Extract the filename from the full path
  filename=$(basename "$file")
  bcftools mpileup -f /nobackup/beegfs/workspace/at991/Data/WI01_2001_043_13_2P.fna "$file" | bcftools call --ploidy 1 -mv -Ob -o calls.bcf
  bcftools view -i '%QUAL>=10' -V indels calls.bcf > "$filename".raw.vcf
  bgzip -f "$filename".raw.vcf
  bcftools index "$filename".raw.vcf.gz
  bcftools consensus -f /nobackup/beegfs/workspace/at991/Data/WI01_2001_043_13_2P.fna "$filename".raw.vcf.gz > "$filename".fasta
done
