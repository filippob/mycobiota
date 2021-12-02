## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
core=8

## Create analysis folders

echo $HOME
echo $currpath

## FastQC
cd $currpath

$HOME/software/FastQC/fastqc /home/ngs/210915_M04028_0141_000000000-JY8M4/*.fastq.gz -o $currpath/Analysis/raw_quality -t 8
cd $currpath/Analysis/raw_quality
multiqc .

echo "DONE!"

