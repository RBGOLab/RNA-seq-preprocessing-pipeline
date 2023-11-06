# Pipeline for aligning pair-ended RNA-seq fastq files to genome

The pipeline consists of a bash script which takes in an specially formatted input file and can automatically perform:

1. Generate a [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) report for the original fastq files
2. Remove adaptor contamination using [Scythe](https://github.com/vsbuffalo/scythe)
3. Remove low quality bases using [Sickle](https://github.com/najoshi/sickle)
4. Align fastq files to the genome using [STAR](https://github.com/alexdobin/STAR)
5. Generate a FASTQC report for the formatted files 
6. Generate a [MULTIQC](https://github.com/ewels/MultiQC) report for formated fastq files 
