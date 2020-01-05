# StrVCTVRE_beta
Early Release version of StrVCTVRE

## To run StrVCTVRE, following these steps:

### 1. Download and install Python (if not done already)
I recommend [installing python through Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

### 2. Download and install the following packages
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

### 3. Download StrVCTVRE
StrVCTVRE requires the 9.2GB [hg38.100way.phyloP100way.bw](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw) file to run. Most users should download StrVCTVRE bundled with hg38.100way.phyloP100way.bw [here](). 

If for some reason you have previously downloaded hg38.100way.phyloP100way.bw and wish to save space, the version of StrVCTVRE without hg38.100way.phyloP100way.bw can be downloaded [here]. You will then need to give StrVCTVRE your the absolute path to your local copy of hg38.100way.phyloP100way.bw as the third argument (see below).

### 4. Unpack StrVCTVRE
