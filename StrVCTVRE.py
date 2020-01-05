
# coding: utf-8

# # This notebook annotates input file SVs using StrVCTVRE
# For a vcf entry to be annotated, it must have an END tag and a SVTYPE tag. Only exonic deletions and duplications will be annotated. Must be in GRCh38. Only annotates autosomes, X, and Y.

# In[55]:

# may need to put each of these in it's own statement, throw an error if fails
import sys
import numpy as np
import pandas as pd
import pybedtools
from cyvcf2 import VCF,Writer
import annotationFinalForStrVCTVRE
reload(annotationFinalForStrVCTVRE)
from joblib import dump, load


# In[53]:

# read in command line arguments
# write appropriate errors also
try:
    vcfPathIn = sys.argv[1]
except:
    sys.exit('Error: no input vcf path supplied')

try:    
    vcfPathOut = sys.argv[2]
except:
    sys.exit('Error: no output vcf path supplied')
try:
    phylopPath = sys.argv[3]
except:
    print('\nusing phyloP file in data subdirectory\n')
    phylopPath = 'data/hg38.phyloP100way.bw'
# vcfPathIn = '/data/andrewsharo-S/thesis/Aim1/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz'
# vcfPathOut = '/data/andrewsharo-S/thesis/Aim1/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.annotatedStrVCTVRE.gz'


# read VCF into one large csv file

# In[48]:

print('\nreading VCF...\n')
toDf = []
toNumpy = []
count = 0
for var in VCF(vcfPathIn,gts012=True):
    if var.INFO.get('END') and var.INFO.get('SVTYPE'):
        count += 1
        if count > 1000:
            break
        entry = np.array([var.CHROM, var.POS, var.INFO['END'], var.INFO['SVTYPE']])
        toDf.append(entry)
df = pd.DataFrame(toDf,columns=['chrom','start','end','svtype'])


# Change formatting, keep only dels and dups

# In[21]:

print('\nformatting VCF data...\n')

# make old index so we can annotate SVs rapidly at the end
df['OldID'] = pd.Series(df.index.values)

# check that the chroms all have chr in front
if sum(df['chrom'].str.startswith('chr',na=False))/df.shape[0] < 0.5:
    df['chrom'] = 'chr' + df['chrom'].astype(str)

acceptedChroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
                 'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
# keep only autosomes, X, and Y 
df = df[df['chrom'].isin(acceptedChroms)].copy()

# keep only deletions and duplications
df = df[((df['svtype'] == 'DEL') | (df['svtype'] == 'DUP'))].copy()
df['DEL'] = df['svtype'] == 'DEL'


# Determine how many exons overlap each variant

# In[23]:

print('\nidentifying exonic deletions and duplications...\n')

exons = pybedtools.BedTool('data/exons_Appris_featurized_transcript_Chr1-Y_loeuf.sorted.bed')
df[['chrom','start','end','OldID']].to_csv('data/svs.bed',sep='\t', index=False,header=False)
a = pybedtools.BedTool('data/svs.bed')
b = a.intersect(exons, wa=True, wb=True).saveas('data/svsExonOverlap.bed')
exonOverlap = pd.read_csv('data/svsExonOverlap.bed', sep='\t', header=None, usecols=[0,1,2,3],
                          names=['chrom', 'start', 'stop', 'OldID'])
exonOverlap['numExons'] = exonOverlap.groupby(by='OldID').chrom.transform('size') # choice of chrom column here is arbitrary
exonOverlap.drop_duplicates(subset='OldID', inplace=True)


# Drop variants that overlap no exons

# In[24]:

out = df.merge(exonOverlap[['numExons','OldID']],how='left',on='OldID')
out = out[out['numExons'] > 0]


# In[26]:

out[['chrom','start','end','OldID','DEL']].to_csv('data/svsForAnnotation.csv')


# Score each variant

# In[32]:

print('\nscoring exonic deletions and duplications...\n')
annotationFinalForStrVCTVRE.annotateSVs('data/svsForAnnotation.csv', 'data/svsAnnotated.csv', phylopPath)


# In[35]:

an = pd.read_csv('data/svsAnnotated.csv')


# In[37]:

rf = load('data/rfMultipleTrainingSets.joblib')

X = an[['DEL','numExonsFinal','phyloP', 'lowestExonRank', 'allSkippable','lowestExonsInGene', 'anyConstExon','pLIMax','loeufMin', 'cdsFracStartMin', 'cdsFracEndMax', 'cdsFracMax', 'pLI_max25_ID', 'loeuf_min25_ID','topExp','topUsage','maxStrength']].copy()

an['path'] = rf.predict_proba(X)[:,1]


# In[56]:

an.set_index('OldID', inplace=True)


# Annotate vcf with StrVCTVRE pathogenicity scores

# In[51]:

print('\nwriting annotated VCF...\n')

vcf = VCF(vcfPathIn)
vcf.add_info_to_header({'ID':'StrVCTVRE','Description':'pathogenicity score for structural variants','Type':'Float','Number':'1'})

w = Writer(vcfPathOut,vcf)

count = 0
for var in vcf:
    if var.INFO.get('END') and var.INFO.get('SVTYPE'):
        if an.index.contains(count):
            var.INFO['StrVCTVRE'] = an.loc[count,'path']
        else:
            var.INFO['StrVCTVRE'] = 'not_exonic'
        count += 1
    else:
        var.INFO['StrVCTVRE'] = 'NA'
    w.write_record(var)
w.close();
vcf.close()


# In[ ]:

print('\nFinished\n')


# In[ ]:



