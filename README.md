## software and steps to analyse mycobiota data

### Software

- fastqc
- multiqc: to combine multiple fastq files (Python code: pip install multiqc)
- cutadapt: to remove barcodes (sequencing primers) [Python code: python3 -m pip install --user --upgrade cutadapt]
- sickle: trimming of reads (C++ code: download and make)
- micca: for i) paired-end reads joining; ii) N-filtering; iii) OTU-picking; iv) assignment of taxonomies; v) phylogenetic tree [Python code: install through docker --> docker pull compmetagen/micca]
- rdp: micca dependency (probably preinstalled in the docker image)

