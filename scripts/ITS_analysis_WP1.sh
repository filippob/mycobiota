

## Set some parameters for the script
currpath=$(pwd) # top project directory
core=8 # number of core to use

## Create analysis folder
mkdir $currpath/raw_quality 
mkdir $currpath/renamed/
mkdir $currpath/cutadapt/
mkdir $currpath/quality_control/
mkdir $currpath/quality_control/cutadapt
mkdir $currpath/quality_control/trimmed
mkdir $currpath/trimmed
mkdir $currpath/MICCA_ITS_WP1


## Visualize quality of raw reads

cd $currpath

fastqc $currpath/Raw_reads/*.fastq.gz -o ./raw_quality -t 8
cd $currpath/raw_quality
multiqc .

## Rename samples


cd $currpath/Raw_reads

este=".fastq.gz"

for i in *.fastq.gz
do
  sample=$(echo "$i" | cut -d "_" -f1 | cut -d "-" -f1,2,3)
  read=$(echo "$i" | cut -d "_" -f4)
  echo $sample
  cp $i $currpath/renamed/$sample"_"$read$este
  echo -e "$i\t-->\t$sample"_"$read$este" >> log_renamer.txt
done

## Count reads

cd $currpath/renamed

for i in *.fastq.gz
do
        echo -n $i >> seq_count_ITS_raw.txt
        echo -n " " >> seq_count_ITS_raw.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_ITS_raw.txt
done

## Create a file of names that will be used for looping. Only file/sample name, remove extension and R1/R2
for i in *.fastq.gz
do 
echo "$i" | cut -d "_" -f1 >> names.txt
sed 'n; d' names.txt > names_single.txt  
done

cp names_single.txt $currpath/trimmed
cp names_single.txt $currpath/cutadapt
cp names_single.txt $currpath/MICCA_ITS_WP1

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
        /usr/local/bin/cutadapt -g Forward=CTTGGTCATTTAGAGGAAGTAA -a RCreverse=GCATCGATGAAGAACGCAGC -o $currpath/cutadapt/"${file}_R1_cutadapt.fastq.gz" "${file}_R1.fastq.gz" 1>> $currpath/quality_control/cutadapt/report_cutadapt_primer_R1.txt
done < names_single.txt


while read file
do
	echo "Running cutadapt on R2 of file "${file}"" 
	/usr/local/bin/cutadapt -g Reverse=GCTGCGTTCTTCATCGATGC -a RCforward=TTACTTCCTCTAAATGACCAAG -o $currpath/cutadapt/"${file}_R2_cutadapt.fastq.gz" "${file}_R2.fastq.gz" 1>> $currpath/quality_control/cutadapt/report_cutadapt_primer_R2.txt
done < names_single.txt



# --discard-untrimmed, --trimmed-only
#                        Discard reads that do not contain an adapter.

#--pair-filter=(any|both|first)
#                        Which of the reads in a paired-end read have to match
#                        the filtering criterion in order for the pair to be
#                        filtered. Default: any

cd $currpath/cutadapt

fastqc ./*.fastq.gz -o $currpath/quality_control/cutadapt -t 8

cd  $currpath/quality_control/cutadapt
multiqc .


# trim low quality part. Q = 20 to be a bit more flexible, ITS has general lower quality than 16S 

cd $currpath/cutadapt

while read file
do
	echo "Running sickle on file "${file}""
	echo "Running sickle on file "${file}"" >> $currpath/quality_control/trimmed/stats_trim.txt
	sickle pe -f "${file}_R1_cutadapt.fastq.gz" -r "${file}_R2_cutadapt.fastq.gz" -o $currpath/trimmed/"${file}_trimmed_R1.fastq.gz" -p $currpath/trimmed/"${file}_trimmed_R2.fastq.gz" -s $currpath/trimmed/"${file}_singles.gz" -t sanger -q 20 -g 1>> $currpath/quality_control/trimmed/stats_trim.txt
done < names_single.txt



cd $currpath/trimmed
fastqc ./*.fastq.gz -o $currpath/quality_control/trimmed -t 8


cd  $currpath/quality_control/trimmed
multiqc .


# Count QC treated data


cd $currpath/trimmed

for i in *.fastq.gz
do
        echo -n $i >> seq_count_ITS_QC.txt
        echo -n " " >> seq_count_ITS_QC.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_ITS_QC.txt
done


## Join reads (MICCA)

cd $currpath/trimmed

# remove singles reads from sickle

rm ./*singles.gz

gunzip *.fastq.gz

micca mergepairs -i $currpath/trimmed/*_R1.fastq -o $currpath/MICCA_ITS_WP1/WP1_assembled_ITS.fastq -l 32 -d 8 -t 7

# -l : minimum overlap between reads 
# -d : maximum mismatch in overlap region

# Counting reads in assembled file

grep -c '^@M' $currpath/MICCA_ITS_WP1/WP1_assembled_ITS.fastq   

# 11534421. Total sum of trimmed reads = 26670364; Teoric 100% assembly = 26670364/2 = 13335182
# Read loss from QC reads to assembled = 13335182 - 11534421 = 1800761;
# Read loss from QC reads to assembled % = 1800761/13335182*100 = 13%

# zipping back trimmed files

cd $currpath/trimmed

gzip *.fastq


# Count assembled reads per sample

cd $currpath/MICCA_ITS_WP1

while read i
do
  echo -n $i >> seq_count_assembled.txt
  echo -n " " >> seq_count_assembled.txt
  grep -wc "$i" WP1_assembled_ITS.fastq >> seq_count_assembled.txt
  echo "##" >> seq_count_assembled.txt
done < names_single.txt

# FASTQC on the assembled to evaluate lenght and quality of assembled
cd $currpath/MICCA_ITS_WP1

fastqc ./WP1_assembled_ITS.fastq -t 8

# Remove N from assembly
micca filter -i $currpath/MICCA_ITS_WP1/WP1_assembled_ITS.fastq -o $currpath/MICCA_ITS_WP1/WP1_assembled_ITS.fasta --maxns 0

# count 
grep -c '>' $currpath/MICCA_ITS_WP1/WP1_assembled_ITS.fasta

# 11440208. Pre N removal was 11534421

## here I should insert some recap file in excel like. All counts are in txt file with space as sep
#raw count: $currpath/renamed/seq_count_16S_raw.txt
#qc count: $currpath/trimmed/seq_count_16S_QC.txt
#assembled: $currpath/MICCA_16S_WP3WP4/seq_count_assembled.txt


# pick otu
micca otu -m denovo_unoise -i $currpath/MICCA_ITS_WP1/WP1_assembled_ITS.fasta -o $currpath/MICCA_ITS_WP1/ -t 8 --rmchim

# classify RDP

export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/
micca classify -m rdp -i $currpath/MICCA_ITS_WP1/otus.fasta --rdp-gene fungalits_unite -o $currpath/MICCA_ITS_WP1/taxa.txt

# also do the phylogenetic tree to be able to use UniFrac or PD measure

micca msa -m muscle -i $currpath/MICCA_ITS_WP1/otus.fasta -o $currpath/MICCA_ITS_WP1/WP1_msa.fasta 

micca tree -i $currpath/MICCA_ITS_WP1/WP1_msa.fasta -o $currpath/MICCA_ITS_WP1/WP1tree.tree -m muscle --muscle-cluster upgma

micca root -i $currpath/MICCA_ITS_WP1/WP1tree.tree -o $currpath/MICCA_ITS_WP1/WP1tree_rooted.tree


