iqtree --prefix Rlowsnps_stripped_nogaps -s Rlow_consensus_snps_trimmed_nogaps.masked.fasta -m TEST -T AUTO -B 1000 --date ~/OneDrive/Data/Cambridge_Project/Edited_Lucy_samples_dates.txt

fasttree -nt -gtr -boot 1000 < Rlow_consensus_snps_trimmed_nogaps.masked.fasta > Fasttree_Rlow_snps_Lucy_trimmed.newick

raxmlHPC-PTHREADS -T 10 -s Rlow_consensus_snps_trimmed_nogaps.masked.fasta -f a -n RaxML_Rlowsnps_stripped_nogaps -m GTRGAMMA4 -c -F -k -b 1234 --bootstop-perms=100