## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
outdir="Analysis/bacteria"
core=8

if [ ! -d "${outdir}/micca_16S" ]; then
	mkdir -p ${outdir}/micca_16S
fi

## OTU picking
# pick otu
singularity run $currpath/micca.sif micca otu -m denovo_unoise -i $currpath/${outdir}/micca_16S/WP1_assembled_16S.fasta -o $currpath/${outdir}/micca_16S/ -t 8 --rmchim

echo "DONE!!"

