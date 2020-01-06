# StrVCTVRE_beta
Early Release version of StrVCTVRE

## To run StrVCTVRE, following these steps:

### 1. Download and install Python (if not done already)
I recommend [installing python through Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

### 2. Download PhyloP conservation scores for human genome 38
This is a 9.2GB file that is required to run StrVCTVRE. It can be downloaded [here](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw). Depending on your connection, it should take about 30 minutes to download.

### 3. Download and install the following packages
* numpy
* pandas
* pybedtools
* cyvcf2
* joblib
* pyBigWig

In conda, numpy can be downloaded and installed using the command 
```
conda install numpy
```
and similarly for the other packages.

### 4. Download StrVCTVRE
Linux users should [download StrVCTVRE tarball](https://github.com/andrewSharo/StrVCTVRE_beta/archive/v.0.1.tar.gz). Other users should download [StrVCTVRE zip file for Windows or Mac](https://github.com/andrewSharo/StrVCTVRE_beta/archive/v.0.1.zip). 

### 5. Uncompress StrVCTVRE
Linux users should run tar -xzf \[filename.tar.gz\]. Windows and Mac users should extract the files from the .zip file. 

Inside the uncompressed folder, you will find several files and a folder called 'data'. Please move the 9.2GB hg38.phyloP100way.bw file to the 'data' folder. If this is not possible, you will need to provide the path to this file when running StrVCTVRE (see below).

### 6. Running StrVCTVRE
To run StrVCTVRE, change directories to the folder containing StrVCTVRE.py, and run 
```
python StrVCTVRE.py path_to_input_vcf_to_annotate path_to_output_annotated_vcf [path_to_hg38.phyloP100way.bw]
``` 
The last argument in optional depending on whether the file is in the 'data' folder. 

For example, a user who moved hg38.phyloP100way.bw to the 'data' folder might run
```
python StrVCTVRE.py /data/patients/patient1.vcf /data/patients/patient1_annotated.vcf
```
A user who has hg38.phyloP100way.bw in a different folder might run
```
python StrVCTVRE.py /data/patients/patient1.vcf /data/patients/patient1_annotated.vcf /data/conservation/hg38.phyloP100way.bw
```

### 7. Other notes
StrVCTVRE can annotate and output vcf files compressed with bgzip. It uses file extensions to know whether an input or output file should be compressed or not.

