#################################################
## Script to rename samples and count input reads
#################################################

## setting the environmnent
currpath=$(pwd)
datapath="/home/ngs/200206_M04028_0114_000000000-CWB58/"
outdir="Analysis/bacteria/renamed"
core=8

## Rename samples
cd $currpath

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

cd $datapath

este=".fastq.gz"

for i in *.fastq.gz
do
  sample=$(echo "$i" | cut -d "_" -f1 | cut -d "-" -f1,2,3)
  read=$(echo "$i" | cut -d "_" -f4)
  echo $sample
  cp $i ${currpath}/${outdir}/$sample"_"$read$este
  echo -e "$i\t-->\t$sample"_"$read$este" >> ${currpath}/${outdir}/log_renamer.txt
done

## Count reads

cd ${currpath}/${outdir}

for i in *.fastq.gz
do
        echo -n $i >> seq_count_16S_raw.txt
        echo -n " " >> seq_count_16S_raw.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_16S_raw.txt
done

echo "DONE!"

