## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
core=8

if [ ! -d "Analysis/micca_its" ]; then
	mkdir Analysis/micca_its
fi


# classify RDP
# export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/
singularity run $currpath/micca.sif micca classify -m rdp -i $currpath/Analysis/micca_its/otus.fasta --rdp-gene 16srrna -o $currpath/Analysis/micca_its/taxa.txt

## CLASSIFY WITH VSEARCH AND SILVA!

# QIIME compatible SILVa DB should be downloaded
#micca classify -m cons -i /$currpath/MICCA_16S_WP1/otus.fasta  -o /$currpath/MICCA_16S_WP1/taxa_SILVA.txt --ref /home/fvitali/Documenti/Silva_132_release/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna --ref-tax /home/fvitali/Documenti/Silva_132_release/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt


echo "DONE!"

