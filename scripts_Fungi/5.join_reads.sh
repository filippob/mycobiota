## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
core=8

if [ ! -d "Analysis/micca_its" ]; then
	mkdir Analysis/micca_its
fi


# Count trimmed data
cd $currpath/Analysis/trimmed

for i in *.fastq.gz
do
        echo -n $i >> seq_count_ITS_QC.txt
        echo -n " " >> seq_count_ITS_QC.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_ITS_QC.txt
done


## Join reads (MICCA)

cd $currpath/Analysis/trimmed

# remove singles reads from sickle

rm ./*singles.gz

gunzip *.fastq.gz

## SINGULARITY CONTAINER ##
singularity run $currpath/micca.sif micca mergepairs -i $currpath/Analysis/trimmed/*_R1.fastq -o $currpath/Analysis/micca_its/WP1_assembled_ITS.fastq -l 32 -d 8 -t 7

# -l : minimum overlap between reads
# -d : maximum mismatch in overlap region

# Counting reads in assembled file

grep -c '^@M' $currpath/Analysis/micca_its/WP1_assembled_ITS.fastq

# 11534421. Total sum of trimmed reads = 26670364; Teoric 100% assembly = 26670364/2 = 13335182
# Read loss from QC reads to assembled = 13335182 - 11534421 = 1800761;
# Read loss from QC reads to assembled % = 1800761/13335182*100 = 13%

# zipping back trimmed files

cd $currpath/Analysis/trimmed

gzip *.fastq

echo "DONE!!"

