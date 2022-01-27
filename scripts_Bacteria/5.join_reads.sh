## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/200206_M04028_0114_000000000-CWB58/"
outdir="Analysis/bacteria"
core=8

if [ ! -d "${outdir}/micca_16S" ]; then
	mkdir -p ${outdir}/micca_16S
fi


# Count trimmed data
cd $currpath/${outdir}/trimmed

echo " - counting sequences "
for i in *.fastq.gz
do
        echo -n $i >> seq_count_16S_QC.txt
        echo -n " " >> seq_count_16S_QC.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_16S_QC.txt
done


## Join reads (MICCA)

cd $currpath/${outdir}/trimmed

# remove singles reads from sickle

if [ -f ./*singles.gz ]; then
	rm ./*singles.gz
fi

echo " - uncompressing trimmed fastq files "
gunzip *.fastq.gz

## SINGULARITY CONTAINER ##
echo " - joining reads"
singularity run $currpath/micca.sif micca mergepairs -i $currpath/${outdir}/trimmed/*_R1.fastq -o $currpath/${outdir}/micca_16S/WP1_assembled_16S.fastq -l 32 -d 8 -t 7

# -l : minimum overlap between reads
# -d : maximum mismatch in overlap region

# Counting reads in assembled file
echo " - counting reads after joining "
grep -c '^@M' $currpath/${outdir}/micca_16S/WP1_assembled_16S.fastq

# 11534421. Total sum of trimmed reads = 26670364; Teoric 100% assembly = 26670364/2 = 13335182
# Read loss from QC reads to assembled = 13335182 - 11534421 = 1800761;
# Read loss from QC reads to assembled % = 1800761/13335182*100 = 13%

# zipping back trimmed files
echo " - recompressing trimmed fastq files "
cd $currpath/${outdir}/trimmed

gzip *.fastq

echo "DONE!!"

