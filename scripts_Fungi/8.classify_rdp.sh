## setting the enviornmnent
currpath=$(pwd)
datapath="/home/ngs/210915_M04028_0141_000000000-JY8M4/"
core=8

if [ ! -d "Analysis/micca_its" ]; then
	mkdir Analysis/micca_its
fi


# classify RDP
# export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/
singularity run $currpath/micca.sif micca classify -m rdp -i $currpath/Analysis/micca_its/otus.fasta --rdp-gene fungalits_unite -o $currpath/Analysis/micca_its/taxa.txt

echo "DONE!"

