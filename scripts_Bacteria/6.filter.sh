## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/200206_M04028_0114_000000000-CWB58/"
outdir="Analysis/bacteria"
core=8

if [ ! -d "${outdir}/micca_16S" ]; then
	mkdir -p ${outdir}/micca_16S
fi

# Remove N from assembly
echo " - filtering reads"
singularity run $currpath/micca.sif micca filter -i $currpath/${outdir}/micca_16S/WP1_assembled_16S.fastq -o $currpath/${outdir}/micca_16S/WP1_assembled_16S.fasta --maxns 0

# count
echo " - countig reads after filtering"
grep -c '>' $currpath/${outdir}/micca_16S/WP1_assembled_16S.fasta

echo "DONE!!"

