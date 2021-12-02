## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
core=8

export PATH=/home/biscarinif/.local/bin:$PATH

## make folders
if [ ! -d "Analysis/cutadapt" ]; then
	mkdir Analysis/cutadapt
fi

if [ ! -d "Analysis/trimmed" ]; then
	mkdir Analysis/trimmed
fi

if [ ! -d "Analysis/micca_its" ]; then
	mkdir Analysis/micca_its
fi

if [ ! -d "Analysis/quality_control" ]; then
	mkdir Analysis/quality_control
fi


## Create a file of names that will be used for looping. Only file/sample name, remove extension and R1/R2

cd $currpath/Analysis/renamed

for i in *.fastq.gz
do
echo "$i" | cut -d "_" -f1 >> names.txt
sed 'n; d' names.txt > names_single.txt
done

cp names_single.txt $currpath/Analysis/trimmed
cp names_single.txt $currpath/Analysis/cutadapt
cp names_single.txt $currpath/Analysis/micca_16S

# remove primers with cutadapt
# Primers (Sequences from Pindo and FMACH)
# forward: CCTACGGGNGGCWGCAG
# reverse: GACTACNVGGGTWTCTAATCC

cd $currpath/renamed

while read file
do
	echo "Running cutadapt on file "${file}""
	/usr/local/bin/cutadapt -g Forward=CCTACGGGNGGCWGCAG -G Reverse=GACTACNVGGGTWTCTAATCC --discard-untrimmed --pair-filter=any -o $currpath/cutadapt/"${file}_R1_cutadapt.fastq.gz" -p $currpath/cutadapt/"${file}_R2_cutadapt.fastq.gz" "${file}_R1.fastq.gz" "${file}_R2.fastq.gz" >> $currpath/quality_control/cutadapt/cutadapt_report.txt  
done < names_single.txt

# --discard-untrimmed, --trimmed-only
#                        Discard reads that do not contain an adapter.

#--pair-filter=(any|both|first)
#                        Which of the reads in a paired-end read have to match
#                        the filtering criterion in order for the pair to be
#                        filtered. Default: any


echo "DONE!!"

