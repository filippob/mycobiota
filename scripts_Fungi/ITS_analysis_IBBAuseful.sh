

## Set some parameters for the script
currpath=$(pwd) # top project directory
core=8 # number of core to use
project=IBBA_useful

## Create analysis folder
mkdir $currpath/raw_quality 
mkdir $currpath/renamed/
mkdir $currpath/cutadapt/
mkdir $currpath/quality_control/
mkdir $currpath/quality_control/cutadapt
mkdir $currpath/quality_control/trimmed
mkdir $currpath/trimmed
mkdir $currpath/MICCA_ITS_$project


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
  sample=$(echo "$i" | cut -d "_" -f1 | cut -d "-" -f1)
  read=$(echo "$i" | cut -d "_" -f4)
  echo "$sample"_"$read$este"
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
cp names_single.txt $currpath/MICCA_ITS_$project

# remove primers with cutadapt

# forward: CTTGGTCATTTAGAGGAAGTAA
# reverse: TCCTCCGCTTATTGATATGC
# RC forward: TTACTTCCTCTAAATGACCAAG
# RC reverse: GCATATCAATAAGCGGAGGA

cd $currpath/renamed

# not using option --discard-untrimmed as not all ITS fragments would have the RCreverse in reads. The option should discard reads without primers, the risk is to loose those

while read file                                                                                        
do
        echo "Running cutadapt on R1 of file "${file}""
        /usr/local/bin/cutadapt -g Forward=CTTGGTCATTTAGAGGAAGTAA -a RCreverse=GCATATCAATAAGCGGAGGA -o $currpath/cutadapt/"${file}_R1_cutadapt.fastq.gz" "${file}_R1.fastq.gz" 1>> $currpath/quality_control/cutadapt/report_cutadapt_primer_R1.txt
done < names_single.txt


while read file
do
	echo "Running cutadapt on R2 of file "${file}"" 
	/usr/local/bin/cutadapt -g Reverse=TCCTCCGCTTATTGATATGC -a RCforward=TTACTTCCTCTAAATGACCAAG -o $currpath/cutadapt/"${file}_R2_cutadapt.fastq.gz" "${file}_R2.fastq.gz" 1>> $currpath/quality_control/cutadapt/report_cutadapt_primer_R2.txt
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

micca mergepairs -i $currpath/trimmed/*_R1.fastq -o $currpath/MICCA_ITS_$project/$project_assembled_ITS.fastq -l 32 -d 8 -t 7

# -l : minimum overlap between reads 
# -d : maximum mismatch in overlap region

# Counting reads in assembled file

grep -c '^@M' $currpath/MICCA_ITS_$project/$project_assembled_ITS.fastq   

# 85003. Total sum of trimmed reads = 19380496; Teoric 100% assembly = 19380496/2 = 9690248
# Read loss from QC reads to assembled = 9690248 - 11534421 = 1800761;
# Read loss from QC reads to assembled % = 1800761/9690248*100 = 13%

# zipping back trimmed files

cd $currpath/trimmed

gzip *.fastq


# Count assembled reads per sample

cd $currpath/MICCA_ITS_$project

while read i
do
  echo -n $i >> seq_count_assembled.txt
  echo -n " " >> seq_count_assembled.txt
  grep -wc "$i" $project_assembled_ITS.fastq >> seq_count_assembled.txt
  echo "##" >> seq_count_assembled.txt
done < names_single.txt

# FASTQC on the assembled to evaluate lenght and quality of assembled
cd $currpath/MICCA_ITS_$project

fastqc ./$project_assembled_ITS.fastq -t 8

# Remove N from assembly
micca filter -i $currpath/MICCA_ITS_$project/$project_assembled_ITS.fastq -o $currpath/MICCA_ITS_$project/$project_assembled_ITS.fasta --maxns 0

# count 
grep -c '>' $currpath/MICCA_ITS_$project/$project_assembled_ITS.fasta

# 11440208. Pre N removal was 11534421

## here I should insert some recap file in excel like. All counts are in txt file with space as sep
#raw count: $currpath/renamed/seq_count_16S_raw.txt
#qc count: $currpath/trimmed/seq_count_16S_QC.txt
#assembled: $currpath/MICCA_16S_WP3WP4/seq_count_assembled.txt



##############################################
######### OTU PICKING ON FORWARD READS  ######

# join reads

cp *_R1.fastq ./R1
cp *_R2.fastq ./R2

cd $currpath/trimmed/R1

micca merge -i $currpath/trimmed/R1/*.fastq -o $currpath/trimmed/R1/merged_R1.fastq

sed -n '1~4s/^@/>/p;2~4p' $currpath/trimmed/R1/merged_R1.fastq >  $currpath/MICCA_ITS_IBBA_useful_R1/IBBA_useful_R1_assembled_ITS.fasta

cd $currpath

# pick otu
micca otu -m denovo_unoise -i $currpath/MICCA_ITS_IBBA_useful_R1/IBBA_useful_R1_assembled_ITS.fasta -o $currpath/MICCA_ITS_IBBA_useful_R1/ -t 8 --rmchim

# classify RDP

export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/
micca classify -m rdp -i $currpath/MICCA_ITS_$project/otus.fasta --rdp-gene fungalits_unite -o $currpath/MICCA_ITS_$project/taxa.txt



##############################################
######### OTU PICKING ON REVERSE READS  ######

# join reads

cd $currpath/trimmed/R2

micca merge -i $currpath/trimmed/R2/*.fastq -o $currpath/trimmed/R2/merged_R2.fastq

sed -n '1~4s/^@/>/p;2~4p' $currpath/trimmed/R2/merged_R2.fastq >  $currpath/MICCA_ITS_IBBA_useful_R2/MICCA_ITS_IBBA_assembled_ITS_R2.fasta


/home/fvitali/Documenti/Personal_PATH_folder/seqtk/seqtk seq -r $currpath/MICCA_ITS_IBBA_useful_R2/MICCA_ITS_IBBA_assembled_ITS_R2.fasta > $currpath/MICCA_ITS_IBBA_useful_R2/MICCA_ITS_IBBA_assembled_ITS_R2_RC.fasta


cd $currpath

# pick otu
micca otu -m denovo_unoise -i $currpath/MICCA_ITS_IBBA_useful_R2/MICCA_ITS_IBBA_assembled_ITS_R2_RC.fasta -o $currpath/MICCA_ITS_IBBA_useful_R2/ -t 5 --rmchim

# classify RDP

export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/
micca classify -m rdp -i $currpath/MICCA_ITS_IBBA_useful_R2/otus.fasta --rdp-gene fungalits_unite -o $currpath/MICCA_ITS_IBBA_useful_R2/taxa.txt


##############################################
######### MAPPING BOTH READS ON THE UNITE ######

# reference

# UNITE https://unite.ut.ee/repository.php (use the general FASTA release)
# Original suggestion/idea  https://www.biostars.org/p/9478987/
# BBmap tool  https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/

# local folder :  /home/fvitali/Documenti/Personal_PATH_folder/BBMap_38.94/bbmap/

# see output of help on the script to use

# 1) index the reference. Should be straightforward, but run times could be long
# 2) map the R1 and R2 on the indexed UNITE. Parameter to look for
#   2a) maxindel=<20>  Don't look for indels longer than this.  Lower is faster.  Set to >=100k for RNA-seq.
#       Maybe in our setting this should be set to 500 or something like that. Consider the result from bioanalyzer on the fragments (see Bianca presentation) we see two peaks with the higher
#       being at approx 850 bp. Maibe a good value could be 1000 bp
#   2b) minratio=<0.56>  minhits=<1> ; setting to higher values should speed up the mapping
# 3) both indexing and mapping should use all avail thread [unless capped with the “t=” flag] 


/home/fvitali/Documenti/Personal_PATH_folder/BBMap_38.94/bbmap/bbsplit.sh build=1 in=IBBA_useful_R1_assembled_ITS.fasta in2=IBBA_useful_R2_assembled_ITS.fasta ref_fungi=sh_general_release_dynamic_10.05.2021.fasta out_fungi=out_fungi.fasta maxindel=1000 t=6  

# maybe is not doing what I want, seems to simply take on the two reda that mapped, but I do not know where


### See if CENTRIFUGE can do; should be similar to Kraken, but more suited for desktop


##############################################
######### Concatenate and pick or blast ######

# This originate from ITSx (https://microbiology.se/software/itsx/). The software removes SSU and 5.8 rRNA genes from a full ITS sequence. The output seems to be a single sequence after 
# extractio  of conserved region, to be used with blast or similar.

# So, the concatenation of R1 and R2 (roughly ITS1 and ITS2 in our case) seems fairly normal and should yeald good results (i.e. the gaps should not introduce enought penalities
# to mess with BLAST identification for example, I am thinking of blast as there is implicit gap dealing, not sure of others like RDP)

# To concatenate, and example is here (article evaluates concatention vs merging in 16S data )

# here the script https://github.com/BEEGlab/16SProcessingPipelinesWithQIIME2/blob/main/TrimPrimersLengthQual_MergeConcat_PE.sh
# he uses PANDAseq to do the concatenation, setting an overlap of 250bp so that most (if not all) the reads fail to be merged, and get concatenated in a single output file. Smart!

sed 's/;.*R1//' merged_R1.fastq > merged_R1_renamed.fastq
sed 's/;.*R2//' merged_R2.fastq > merged_R2_renamed.fastq


pandaseq -f merged_R1_renamed.fastq\
         -r merged_R2_renamed.fastq\
         -o 200\
         -g pandaseq_concatenation.log.txt\
         -w output_concat.fastq\
         -U output_unaligned.fastq -B -T 4

# ok, file "output_unalignes.fasta"

# Before BLAST, should select only the sequences that were used to build an OTU. Can use grep, one liner solution did not work, using a while read do loop


while read seq
do
	grep -A 1 $seq output_unaligned.fasta >> grep_selected_seqs.fasta
done < outids.txt

 grep ">" grep_selected_seqs.fasta | wc -l

 # found only 66 sequence in the concatenated file (out of 246 OTUs): PANDASEQ has lost some sequences! Why?? maybe some default param of quality or lenght or similar
 # should try to do it manually with grep and while read do, so I am sure

# convert to fasta to have less lines
sed -n '1~4s/^@/>/p;2~4p' merged_R1_renamed.fastq >  merged_R1_renamed.fasta
sed -n '1~4s/^@/>/p;2~4p' merged_R2_renamed.fastq >  merged_R2_renamed.fasta

# obtain RC of R2
/home/fvitali/Documenti/Personal_PATH_folder/seqtk/seqtk seq -r merged_R2_renamed.fasta > merged_R2_renamed_RC.fasta

#need to remove las :*** from otuids, number 8

while read ids
do 
echo "$ids" | cut -d ":" -f1,2,3,4,5,6,7 >> otuids_renamed.txt
done < otuids.txt

# loop to concatenate
while read seq
do
	R1=$(grep -A 1 "$seq" merged_R1_renamed.fasta | sed "1 d")
	#echo $R1
	R2=$(grep -A 1 "$seq" merged_R2_renamed_RC.fasta | sed "1 d") 
	#echo $R2
	concat="$R1"-"$R2" 
	echo ">""$seq" >> output.fasta
	echo $concat >> output.fasta
	#write 
done < otuids_renamed.txt

grep ">" output.fasta | wc -l
#246 !
