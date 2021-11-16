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
cp names_single.txt $currpath/Analysis/micca_its

# remove primers with cutadapt
# Primers (Sequences from Pindo and FMACH) consider that in ITS we always remove F and R primer from both reads

# forward: CTTGGTCATTTAGAGGAAGTAA
# reverse: GCTGCGTTCTTCATCGATGC
# RC forward: TTACTTCCTCTAAATGACCAAG
# RC reverse: GCATCGATGAAGAACGCAGC

cd $currpath/renamed

# not using option --discard-untrimmed as not all ITS fragments would have the RCreverse in reads. The option should discard reads without primers, the risk is to loose those

while read file
do
        echo "Running cutadapt on R1 of file "${file}""
        cutadapt -g Forward=CTTGGTCATTTAGAGGAAGTAA -a RCreverse=GCATCGATGAAGAACGCAGC -o $currpath/Analysis/cutadapt/"${file}_R1_cutadapt.fastq.gz" "${file}_R1.fastq.gz" 1>> $currpath/Analysis/quality_control/report_cutadapt_primer_R1.txt
done < names_single.txt


while read file
do
	echo "Running cutadapt on R2 of file "${file}""
	cutadapt -g Reverse=TCCTCCGCTTATTGATATGC -a RCforward=GCATATCAATAAGCGGAGGA -o $currpath/Analysis/cutadapt/"${file}_R2_cutadapt.fastq.gz" "${file}_R2.fastq.gz" 1>> $currpath/Analysis/quality_control/report_cutadapt_primer_R2.txt
done < names_single.txt

echo "DONE!!"

