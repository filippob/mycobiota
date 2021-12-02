## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
core=8

if [ ! -d "Analysis/micca_its" ]; then
	mkdir Analysis/micca_its
fi

# Remove N from assembly
singularity run $currpath/micca.sif micca filter -i $currpath/Analysis/micca_its/WP1_assembled_ITS.fastq -o $currpath/Analysis/micca_its/WP1_assembled_ITS.fasta --maxns 0

# count
grep -c '>' $currpath/Analysis/micca_its/WP1_assembled_ITS.fasta

echo "DONE!!"

