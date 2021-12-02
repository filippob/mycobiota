## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
core=8

if [ ! -d "Analysis/micca_its" ]; then
	mkdir Analysis/micca_its
fi

## OTU picking
# pick otu
singularity run $currpath/micca.sif micca otu -m denovo_unoise -i $currpath/Analysis/micca_its/WP1_assembled_ITS.fasta -o $currpath/Analysis/micca_its/ -t 8 --rmchim

echo "DONE!!"

