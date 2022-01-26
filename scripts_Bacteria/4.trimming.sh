## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/200206_M04028_0114_000000000-CWB58/"
outdir="Analysis/bacteria"
core=8
sickle_exe="$HOME/software/sickle/sickle"

if [ ! -d "${outdir}/trimmed" ]; then
	mkdir -p ${outdir}/trimmed
fi

## TRIMMING
# trim low quality part. 
# Q = 20 inspect quality, eventually for 16S Q can be set to 25

cd ${currpath}/${outdir}/cutadapt

while read file
do
	echo "Running sickle on file "${file}""
	echo "Running sickle on file "${file}"" >> ${currpath}/${outdir}/quality_control/stats_trim.txt
	$sickle_exe pe -f "${file}_R1_cutadapt.fastq.gz" -r "${file}_R2_cutadapt.fastq.gz" -o ${currpath}/${outdir}/trimmed/"${file}_trimmed_R1.fastq.gz" -p ${currpath}/${outdir}/trimmed/"${file}_trimmed_R2.fastq.gz" -s ${currpath}/${outdir}/trimmed/"${file}_singles.gz" -t sanger -q 20 -g 1>> ${currpath}/${outdir}/quality_control/stats_trim.txt
done < names_single.txt

echo "DONE!!"
