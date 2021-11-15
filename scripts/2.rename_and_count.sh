## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
core=8

## Rename samples
cd $currpath

if [ ! -d "Analysis/renamed" ]; then
	mkdir Analysis/renamed
fi

cd $datapath

este=".fastq.gz"

for i in *.fastq.gz
do
  sample=$(echo "$i" | cut -d "_" -f1 | cut -d "-" -f1,2,3)
  read=$(echo "$i" | cut -d "_" -f4)
  echo $sample
  cp $i $currpath/Analysis/renamed/$sample"_"$read$este
  echo -e "$i\t-->\t$sample"_"$read$este" >> $currpath/Analysis/renamed/log_renamer.txt
done

## Count reads

cd $currpath/Analysis/renamed

for i in *.fastq.gz
do
        echo -n $i >> seq_count_ITS_raw.txt
        echo -n " " >> seq_count_ITS_raw.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_ITS_raw.txt
done

echo "DONE!"
