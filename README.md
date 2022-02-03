![StrVCTVRE logo](images/StrVCTVRELogoRandom.PNG)

Structural variant impact predictor developed by Andrew Sharo, Zhiqiang Hu, and Steven Brenner at UC Berkeley and Shamil R. Sunyaev at Harvard Medical School. StrVCTVRE annotates exonic deletions and duplications with a score from 0 to 1, where higher scores are more likely to be pathogenic. Accepts vcf and bed file input in both GRCh38 and GRCh37. StrVCTVRE learns to identify pathogenic variants using a random forest framework trained on data from ClinVar, gnomAD, and a recent great ape sequencing study. StrVCTVRE stands for <ins>Str</ins>uctural <ins>V</ins>ariant <ins>C</ins>lassifier <ins>T</ins>rained on <ins>V</ins>ariants <ins>R</ins>are and <ins>E</ins>xonic.

### StrVCTVRE is now published in the [American Journal of Human Genetics](https://doi.org/10.1016/j.ajhg.2021.12.007) 

### \*\*New\*\*: Visit our [Web Server](https://strvctvre.berkeley.edu) to annotate vcf or bed files with StrVCTVRE scores, or query a single SV

## To run StrVCTVRE, follow these steps:

### 1. Download and install Python (if not done already)
I recommend [installing python through Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

### 2. Clone github StrVCTVRE repository
If [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) is installed on your system, run
```
git clone https://github.com/andrewSharo/StrVCTVRE
```
which will copy the StrVCTVRE files onto your system. Continue to step 3.

If you cannot clone the github repository, Linux users should download the [StrVCTVRE tarball](https://github.com/andrewSharo/StrVCTVRE/archive/v.1.7.tar.gz). Other users should download the [StrVCTVRE zip file for Windows or Mac](https://github.com/andrewSharo/StrVCTVRE/archive/v.1.7.zip). 

To extract the files, linux users should run tar -xzf \[filename.tar.gz\]. Windows and Mac users should extract the files from the .zip file. 

### 3. Download PhyloP conservation scores for human genome 38
This is a 9.2GB file that is required to run StrVCTVRE. It can be downloaded [here](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw). Depending on your connection, it should take about 20 minutes to download.

Please move this 9.2GB hg38.phyloP100way.bw file to the 'data' folder inside the cloned StrVCTVRE repository. If this is not possible, you will need to provide the path to this file when running StrVCTVRE (see below).

### 4. Download and install required packages
If using conda, navigate to the StrVCTVRE folder and run
```
conda env create -f environment_py2.7.yml
```
which will install all necessary packages in a new conda environment. Next run
```
conda activate StrVCTVRE_py_2.7
```
which will activate this new conda environment. If this is successful, continue to Step 5.

If you were not able to use the yml file to create a new conda environment, you will need to manually install the following packages:

* numpy
* pandas
* joblib
* scikit-learn
* pybedtools
* cyvcf2
* pybigwig

In conda, these packages can be downloaded and installed using the command 
```
conda install numpy pandas joblib scikit-learn
```
Some packages are distributed through bioconda, so they can be downloaded and install using the command
```
conda install -c bioconda pybedtools cyvcf2 pybigwig
```

### 5. Test StrVCTVRE for GRCh38
To test that StrVCTVRE is annotating correctly, change your current working directory to the StrVCTVRE folder that contains test_StrVCTVRE.py and run 
```
python test_StrVCTVRE.py 
```
Note that you will need to use the -p flag to provide the path to hg38.phyloP100way.bw if it is not in 'data' folder (see details below). For example, you might run
```
python test_StrVCTVRE.py -p /home/conservation/hg38.phyloP100way.bw
```
If this function prints "SUCCESS" then your copy of StrVCTVRE is working correctly. You may disregard warnings. If the function prints "ERROR" then please [raise a new issue](https://github.com/andrewSharo/StrVCTVRE/issues)

### 6. Test StrVCTVRE for GRCh37
You may skip this section if your structural variants are already in GRCh38.

StrVCTVRE can annotate structural variants in GRCh37. This requires having a functioning copy of the LiftOver software on your system. LiftOver can be downloaded from [the UCSC Genome Browser store](https://genome-store.ucsc.edu/) for free for academic researchers. You do not need to download any liftOver chain.gz files, as they are already included in the StrVCTVRE repository. Once you have downloaded LiftOver onto your system, test StrVCTVRE by running
```
python test_StrVCTVRE_GRCh37.py -l /path/to/liftOver 
```
If this function prints "SUCCESS" then your copy of StrVCTVRE is ready to annotate GRCh37 structural variants. You may disregard warnings. If the function prints "ERROR" then please [raise a new issue](https://github.com/andrewSharo/StrVCTVRE/issues)

### 7. Run StrVCTVRE for GRCh38
To run StrVCTVRE, change your working directory to the folder containing StrVCTVRE.py, and run 
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
python StrVCTVRE.py -i /home/user/patient1.vcf -o /home/user/patient1_annotated.vcf -p /home/conservation/hg38.phyloP100way.bw
```
A user who wants to annotate a bed file and who moved hg38.phyloP100way.bw to the 'data' folder might run
```
python StrVCTVRE.py -i /home/user/patient1.bed -o /home/user/patient1_annotated.bed -f bed 
```

### 8. Run StrVCTVRE for GRCh37
The commands to run StrVCTVRE for GRCh37 are identical to the above commands for GRCh38, except you must explicitly include an argument indicating the assembly is GRCh37, and you must also include the path to the LiftOver executable. For example,
```
python StrVCTVRE.py -i /home/user/patient1.vcf -o /home/user/patient1_annotated.vcf -a GRCh37 -l /path/to/LiftOver
```
Although the SVs are in GRCh37, you will still need to download the phyloP hg38.phyloP100way.bw file. Unless you have put that file in the 'data' folder, you will need to indicate the path to that file, as done here
```
python StrVCTVRE.py -i /home/user/patient1.vcf -o /home/user/patient1_annotated.vcf -a GRCh37 -l /path/to/LiftOver -p /home/conservation/hg38.phyloP100way.bw
```

### 9. Other notes
StrVCTVRE can annotate both vcf and bed files. For a vcf record to be annotated by StrVCTVRE, the record must overlap an exon, have 'END' and 'SVTYPE' entries in the INFO column, and SVTYPE must be 'DUP' or 'DEL'. SVs that do not meet these requirements will be annotated with a string describing what they are missing. Deletions and Duplications larger than 3Mb will be given a score of 1.0 since very few benign SVs are larger than 3Mb. StrVCTVRE can annotate and output vcf files compressed with bgzip. It uses file extensions to know whether an input or output file should be compressed or not. 

For a bed file to be annotated by StrVCTVRE it must have no header and match this format:

chromosome[tab]start[tab]end[tab]svtype

See below for an example file. 
```
chr1  123456  234567  DEL
chr2  345678  456789  DUP
```
If you have any questions that aren't answered here, please [raise a new issue](https://github.com/andrewSharo/StrVCTVRE/issues)

## Citation
If you use StrVCTRE in your work, please cite:

Sharo AG, Hu Z, Sunyaev SR, Brenner SE. 2022. StrVCTVRE: A supervised learning method to predict the pathogenicity of human genome structural variants. *The American Journal of Human Genetics 109*:195-209. doi:[10.1016/j.ajhg.2021.12.007](https://doi.org/10.1016/j.ajhg.2021.12.007)

