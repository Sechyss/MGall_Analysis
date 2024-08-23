#!/bin/bash -l

#SBATCH --partition=defq
#SBATCH --job-name=VCF_Rlow
#SBATCH --time=96:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks=32                         # Run a single task (increase this value for parallelization across CPUs)
#SBATCH --account=c.bonneaud                # The accounting code - usually named after the PI for a project
#SBATCH --constraint=IB                     # Specify features required by the job
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

cd /nobackup/beegfs/workspace/at991/Data/Mapped_output_Rlow/ || exit

for file in /nobackup/beegfs/workspace/at991/Data/Mapped_output_Rlow/*bam; do
  # Extract the filename from the full path
  filename=$(basename "$file")
  if [ -e "$filename".bam.raw.vcf ]; then
        echo "Output file $filename already exists. Skipping $filename."
        continue
  fi
  bcftools mpileup -f /nobackup/beegfs/workspace/at991/Data/R.fna --threads 32 "$file" | bcftools call --threads 32 --ploidy 1 -mv -Ob -o calls.bcf
  bcftools view --threads 40 -i '%QUAL>=10' -V indels calls.bcf > "$filename".raw.vcf
  bgzip -f "$filename".raw.vcf
  bcftools index --threads 20 "$filename".raw.vcf.gz
  bcftools consensus -f /nobackup/beegfs/workspace/at991/Data/R.fna "$filename".raw.vcf.gz > "$filename".fasta
done
