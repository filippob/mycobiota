## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
core=8
sickle_exe="$HOME/software/sickle/sickle"

if [ ! -d "Analysis/trimmed" ]; then
	mkdir Analysis/trimmed
fi

## TRIMMING
# trim low quality part. Q = 20 to be a bit more flexible, ITS has general lower quality than 16S

cd $currpath/Analysis/cutadapt

while read file
do
	echo "Running sickle on file "${file}""
	echo "Running sickle on file "${file}"" >> $currpath/Analysis/quality_control/stats_trim.txt
	$sickle_exe pe -f "${file}_R1_cutadapt.fastq.gz" -r "${file}_R2_cutadapt.fastq.gz" -o $currpath/Analysis/trimmed/"${file}_trimmed_R1.fastq.gz" -p $currpath/Analysis/trimmed/"${file}_trimmed_R2.fastq.gz" -s $currpath/Analysis/trimmed/"${file}_singles.gz" -t sanger -q 20 -g 1>> $currpath/Analysis/quality_control/stats_trim.txt
done < names_single.txt

echo "DONE!!"
