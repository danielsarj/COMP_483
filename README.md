# COMP 383/483 - Python wrapper mini-project

## Overview
This repo contains a python wrapper developed for the COMP 383/483 class of 2022/1. Briefly, this script retrieves [SRX5005282](https://www.ncbi.nlm.nih.gov/sra/SRX5005282) Illumina reads; then, using the aforementioned reads, a genome is assembled and annotated. Additionally, using RNA-Seq data from the [SRR1411276](https://www.ncbi.nlm.nih.gov/sra/SRX604287) project and the gene annotation data, gene expression is quantified.  

## Dependencies
It's recommended that the user adds all dependencies to their path or place their respective folders in the working directory. Else, the script might not work properly. 

* [Biopython](https://biopython.org/)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
* [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html)
* [GeneMarkS-2](http://exon.gatech.edu/GeneMark/index.html)
* [Samtools](http://www.htslib.org/)
* [SPAdes](https://github.com/ablab/spades)
* [SRA-Toolkit](https://github.com/ncbi/sra-tools)
* [TopHat2](https://ccb.jhu.edu/software/tophat/index.shtml)


## Instructions
Copy the `miniproject_wrap.py` and the `Ecoli.fasta` files in this repo to the folder which will be the working directory. If dependencies were not added to the user's path, place their folders in the working directory as well. Then, to run the script, use the command:
```
python3 miniproject_wrap.py
```
All outputs generated will be written into a subfolder named `results`. This subfolder will be inside of the working directory. 

## Credits
Python script written by [Daniel S. Araújo](https://orcid.org/0000-0001-6219-9463).
