# StrVCTVRE
Structural variant impact predictor developed by Andrew Sharo and Steven Brenner at UC Berkeley. StrVCTVRE annotates exonic deletions and duplications with a score from 0 to 1, where higher scores are more likely to be pathogenic. Accepts vcf and bed file input. Manuscript in development. StrVCTVRE stands for <ins>Str</ins>uctural <ins>V</ins>ariant <ins>C</ins>lassifier <ins>T</ins>rained on <ins>V</ins>ariants <ins>R</ins>are and <ins>E</ins>xonic.

## To run StrVCTVRE, follow these steps:

### 1. Download and install Python (if not done already)
I recommend [installing python through Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

### 2. Download PhyloP conservation scores for human genome 38
This is a 9.2GB file that is required to run StrVCTVRE. It can be downloaded [here](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw). Depending on your connection, it should take about 30 minutes to download.

### 3. Download and install the following packages
* numpy
* pandas
* joblib
* scikit-learn
* pybedtools
* cyvcf2
* pybigwig

In conda, numpy can be downloaded and installed using the command 
```
conda install numpy
```
and similarly for pandas, joblib, and scikit-learn.

pybedtools is distributed through bioconda, so it can be downloaded and install using the command
```
conda install -c bioconda pybedtools
```
and similarly for cyvcf and pybigwig.

### 4. Download StrVCTVRE
Linux users should download the [StrVCTVRE tarball](https://github.com/andrewSharo/StrVCTVRE/archive/v.1.1.tar.gz). Other users should download the [StrVCTVRE zip file for Windows or Mac](https://github.com/andrewSharo/StrVCTVRE/archive/v.1.1.zip). 

### 5. Uncompress StrVCTVRE
Linux users should run tar -xzf \[filename.tar.gz\]. Windows and Mac users should extract the files from the .zip file. 

Inside the uncompressed folder, you will find several files and a folder called 'data'. Please move the 9.2GB hg38.phyloP100way.bw file to the 'data' folder. If this is not possible, you will need to provide the path to this file when running StrVCTVRE (see below).

### 6. Test StrVCTVRE
To test that StrVCTVRE is annotating correctly, change your current working directory to the uncompressed folder that contains test_StrVCTVRE.py and run 
```
python test_StrVCTVRE.py 
```
Note that you will need to use the -p flag to provide the path to hg38.phyloP100way.bw if it is not in 'data' folder (see details below). For example, you might run
```
python test_StrVCTVRE.py -p /home/conservation/hg38.phyloP100way.bw
```
If this function prints "SUCCESS" then your copy of StrVCTVRE is working correctly. You may disregard warnings. If the function prints "ERROR" then please [raise a new issue](https://github.com/andrewSharo/StrVCTVRE/issues)
### 7. Run StrVCTVRE
To run StrVCTVRE, change your working directory to the uncompressed folder containing StrVCTVRE.py, and run 
```
python StrVCTVRE.py -i /path/to/input/file -o /path/to/output/file [-f {vcf,bed}] [-p path/to/hg38.phyloP100way.bw]
``` 
The -f argument is optional and defaults to vcf. The -p argument is optional depending on whether hg38.phyloP100way.bw is in the 'data' folder. 

For example, a user who wants to annotate a vcf file and who moved hg38.phyloP100way.bw to the 'data' folder might run
```
python StrVCTVRE.py -i /home/user/patient1.vcf -o /home/user/patient1_annotated.vcf 
```
A user who has hg38.phyloP100way.bw in a different folder might run
```
python StrVCTVRE.py /home/user/patient1.vcf /home/user/patient1_annotated.vcf /home/conservation/hg38.phyloP100way.bw
```
A user who wants to annotate a bed file and who moved hg38.phyloP100way.bw to the 'data' folder might run
```
python StrVCTVRE.py -i /home/user/patient1.bed -o /home/user/patient1_annotated.bed -f bed 
```

### 8. Other notes
StrVCTVRE can annotate both vcf and bed files. For a vcf record to be annotated by StrVCTVRE, the record must overlap an exon, have 'END' and 'SVTYPE' entries in the INFO column, and SVTYPE must be 'DUP' or 'DEL'. SVs that do not meet these requirements will be annotated with a string describing what they are missing. StrVCTVRE can annotate and output vcf files compressed with bgzip. It uses file extensions to know whether an input or output file should be compressed or not. 

For a bed file to be annotated by StrVCTVRE it must have no header and match this format:chromosome[tab]start[tab]end[tab]svtype. See below for an example file. 
```
chr1  123456  234567  DEL
chr2  345678  456789  DUP
```
If you have any questions that aren't answered here, please submit a new issue.
