## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/200206_M04028_0114_000000000-CWB58/"
outdir="Analysis/bacteria"
core=8
## primer sequences (primers follow the adapters)
fwd_primer="CCTACGGGNGGCWGCAG" #(adapter: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG)
rev_primer="GACTACHVGGGTATCTAATCC" #(adapter: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG)

export PATH=/home/biscarinif/.local/bin:$PATH

## make folders
if [ ! -d "${outdir}/cutadapt" ]; then
	mkdir -p ${outdir}/cutadapt
fi

if [ ! -d "${outdir}/trimmed" ]; then
	mkdir -p ${outdir}/trimmed
fi

if [ ! -d "${outdir}/micca_its" ]; then
	mkdir -p ${outdir}/micca_its
fi

if [ ! -d "${outdir}/quality_control" ]; then
	mkdir -p ${outdir}/quality_control
fi


## Create a file of names that will be used for looping. Only file/sample name, remove extension and R1/R2

cd "${currpath}/${outdir}/renamed"

for i in *.fastq.gz
do
echo "$i" | cut -d "_" -f1 >> names.txt
sed 'n; d' names.txt > names_single.txt
done

cp names_single.txt ${currpath}/${outdir}/trimmed
cp names_single.txt ${currpath}/${outdir}/cutadapt
cp names_single.txt ${currpath}/${outdir}/micca_16S

# remove primers with cutadapt
# Primers (Sequences from Pindo and FMACH)
# forward: CCTACGGGNGGCWGCAG
# reverse: GACTACNVGGGTWTCTAATCC

cd "${currpath}/${outdir}/renamed"
echo $(pwd)

while read file
do
	echo "Running cutadapt on file '${file}'"
	cutadapt -g Forward=${fwd_primer} -G Reverse=${rev_primer} --discard-untrimmed --pair-filter=any -o "${currpath}/${outdir}/cutadapt/${file}_R1_cutadapt.fastq.gz" -p "${currpath}/${outdir}/cutadapt/${file}_R2_cutadapt.fastq.gz" "${file}_R1.fastq.gz" "${file}_R2.fastq.gz" >> "${currpath}/${outdir}/quality_control/cutadapt_report.txt"  

done < names_single.txt

# --discard-untrimmed, --trimmed-only
#                        Discard reads that do not contain an adapter.

#--pair-filter=(any|both|first)
#                        Which of the reads in a paired-end read have to match
#                        the filtering criterion in order for the pair to be
#                        filtered. Default: any


echo "DONE!!"

