## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/200206_M04028_0114_000000000-CWB58/"
outdir="Analysis/bacteria/raw_quality"
core=8

## Create analysis folders

echo $HOME
echo $currpath

## FastQC
cd $currpath

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

$HOME/software/FastQC/fastqc ${datapath}*.fastq.gz -o ${currpath}/${outdir} -t 8
cd ${currpath}/${outdir}
multiqc .

echo "DONE!"

