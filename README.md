## software and steps to analyse mycobiota data

We are referring here to paired-end reads from ITS sequencing

### Software

- fastqc
- multiqc: to combine multiple fastq files (Python code: pip install multiqc)
- cutadapt: to remove barcodes (sequencing primers) [Python code: python3 -m pip install --user --upgrade cutadapt]
- sickle: trimming of reads (C++ code: download and make)
- micca: for i) paired-end reads joining; ii) N-filtering; iii) OTU-picking; iv) assignment of taxonomies; v) phylogenetic tree [Python code: install through docker --> docker pull compmetagen/micca]
- rdp: micca dependency (probably preinstalled in the docker image)

### Steps

- activate the conda env: `conda activate mycobiota`
- fastqc --> multiqc (`1.fastqc.sh`)
- rename and count reads (`2.rename_and_count.sh`)
- cutadapt (WARNING: The script cutadapt is installed in '/home/biscarinif/.local/bin' which is not on PATH. Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.) (`3.remove_adapters.sh`)
- trim reads (`4.trimming.sh`)
- joining paired-end reads (`5.join_reads.sh`)
- filtering (removing N's) (`6.filter.sh`)
- OTU picking (binning quasi-identical bins) (`7.otu_picking.sh`)
- RDP classifier (assign OTUs to taxonomical classification) (`8.classify_rdp.sh`)

