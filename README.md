# SysISPD
SysFate Illumina Spatial Promoters Demultiplexer

## Parameters

- --gtf : GTF file to use with feature count

- --pos : TSV file who contain coordinates code with their position in the grid

- --samdir : Directory of SAM files

### Optional ones : 

- --modgtf n : n number of kilobases from the promoters

- --merge : to merge all sam files

- --quantifie : parameters for feature count - by default : "featureCounts -T __core__ -F GTF -t __feature__ -s 0 -a __GTF__ -o __out__ -O"

- --feature : features used with feature count - by default : "5UTR,transcript"

## Feature Count Explanation

Feature Count is a read summarization program suitable for counting reads generated from either RNA or genomic DNA sequencing experiments. It works with either single or paired-end reads and provides a wide range of options appropriate for different sequencing applications. FeatureCounts takes as input SAM/BAM files and an annotation file (GTF) including chromosomal coordinates of features. It outputs numbers of reads assigned to features (or meta-features). It also outputs stat info for the overall summrization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons (these reasons are included in the stat info). 

## Script Explanation

The script launches Feature Count on the .sam files in the given folder, adjusting an interval of size nkb around the promoters of the sequences provided by the gtf file, n being given by the --modgtf parameter. The results of the feature count are stored in a file named after the execution date (modifiable with the --name argument) in the current folder (modifiable with --out).
