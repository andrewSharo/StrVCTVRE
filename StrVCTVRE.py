
# coding: utf-8

# # This notebook annotates input file SVs using StrVCTVRE
# For a vcf entry to be annotated, it must have an END tag and a SVTYPE tag. Only exonic deletions and duplications will be annotated. Must be in GRCh38. Only annotates autosomes, X, and Y.

# In[23]:

# may need to put each of these in it's own statement, throw an error if fails
import sys
import numpy as np
import pandas as pd
import pybedtools
from cyvcf2 import VCF,Writer
import annotationFinalForStrVCTVRE
from joblib import dump, load
import argparse
import tempfile
import shutil
import os


# In[24]:

parser = argparse.ArgumentParser(description='Annotate the pathogenicity of exonic deletions and duplications in GRCh38.')
parser.add_argument('-i','--input',help='Input file path',required=True,metavar = '/path/to/input/file',dest='pathIn')
parser.add_argument('-o','--output',help='Output file path',required=True,metavar = '/path/to/output/file',dest='pathOut')
parser.add_argument('-f','--format',help='Input file format, either vcf or bed, defaults to vcf when not provided',choices=['vcf','bed'],dest='formatIn',default='vcf')
parser.add_argument('-p','--phyloP',help='phyloP file path, defaults to \'data/hg38.phyloP100way.bw\' when not provided',default='data/hg38.phyloP100way.bw',
                    metavar = 'path/to/hg38.phyloP100way.bw',dest='phylopPath')
#parser.add_argument('-a','--assembly',help='Genome assembly, either GRCh38 or GRCh37',choices=['GRCh37','GRCh38'])
args = parser.parse_args()


# Create temporary directory to store files created, deleted after finished running

# In[25]:

td = tempfile.mkdtemp(prefix='StrVCTVRE.',suffix='.tmp')


# read VCF or BED into one large csv file

# In[26]:

# if VCF
if args.formatIn == 'vcf':
    print('\nreading VCF...\n')
    toDf = []
    for var in VCF(args.pathIn,gts012=True):
        if var.INFO.get('END') and var.INFO.get('SVTYPE'):
            entry = np.array([var.CHROM, var.POS, var.INFO['END'], var.INFO['SVTYPE']])
            toDf.append(entry)
    df = pd.DataFrame(toDf,columns=['chrom','start','end','svtype'])

#if BED    
else:
    print('\nreading BED...\n')
    toDf = []
    df = pd.read_csv(args.pathIn,sep='\t',names=['chrom','start','end','svtype'],header=None)
    


# Check bed file input has SVTYPE, an easy thing to forget

# In[ ]:

if df['svtype'].isnull().all() & (args.formatIn == 'bed'):
    sys.exit('ERROR: likely missing SVTYPE column from bed file')


# Change formatting, keep only dels and dups

# In[27]:

print('\nformatting VCF data...\n')

# make old index so we can annotate SVs rapidly at the end
df['OldID'] = pd.Series(df.index.values)

# check that the chroms all have chr in front
if sum(df['chrom'].astype(str).str.startswith('chr',na=False))/df.shape[0] < 0.5:
    df['chrom'] = 'chr' + df['chrom'].astype(str)

acceptedChroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
                 'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
# keep only autosomes, X, and Y 
df = df[df['chrom'].isin(acceptedChroms)].copy()
validChrom = df.copy()
validChrom['validChrom'] = True

# remove all start and end values that are not numeric
df = df[pd.to_numeric(df['start'], errors='coerce').notnull()].copy()
df = df[pd.to_numeric(df['end'], errors='coerce').notnull()].copy()
# convert from string to float (relevant to vcf only)
df['start'] = df['start'].astype(float)
df['end'] = df['end'].astype(float)
# check all start and end values are integers
df = df[df['start'] == df['start'] // 1]
df = df[df['end'] == df['end'] // 1]
df['start'] = df['start'].astype(int)
df['end'] = df['end'].astype(int)
validStartEnd = df.copy()

# keep only SVs 50bp or longer
df['length'] = df['end'].astype(int) - df['start'].astype(int)
df = df[df['length'] > 49].copy()
validLength = df.copy()
validLength['validLength'] = True

# keep only deletions and duplications
df = df[((df['svtype'] == 'DEL') | (df['svtype'] == 'DUP'))].copy()
validSVType = df.copy()
validSVType['validSVType'] = True
df['DEL'] = df['svtype'] == 'DEL'


# Determine how many exons overlap each variant

# In[29]:

print('\nidentifying exonic deletions and duplications...\n')

exons = pybedtools.BedTool('data/exons_Appris_featurized_transcript_Chr1-Y_loeuf.sorted.bed')
df[['chrom','start','end','OldID']].to_csv(os.path.join(td,'svs.bed'),sep='\t', index=False,header=False)
a = pybedtools.BedTool(os.path.join(td,'svs.bed'))
b = a.intersect(exons, wa=True, wb=True).saveas(os.path.join(td,'svsExonOverlap.bed'))
exonOverlap = pd.read_csv(os.path.join(td,'svsExonOverlap.bed'), sep='\t', header=None, usecols=[0,1,2,3],
                          names=['chrom', 'start', 'stop', 'OldID'])
exonOverlap['numExons'] = exonOverlap.groupby(by='OldID').chrom.transform('size') # choice of chrom column here is arbitrary
exonOverlap.drop_duplicates(subset='OldID', inplace=True)


# Drop variants that overlap no exons

# In[30]:

out = df.merge(exonOverlap[['numExons','OldID']],how='left',on='OldID')
out = out[out['numExons'] > 0]
validExon = out.set_index('OldID').copy()
validExon['validExon'] = True
# only annotate vars less than 3Mb
out = out[out['length'] < 3000000]


# In[32]:

out[['chrom','start','end','OldID','DEL']].to_csv(os.path.join(td,'svsForAnnotation.csv'))


# Score each variant

# In[33]:

print('\nscoring exonic deletions and duplications...\n')
annotationFinalForStrVCTVRE.annotateSVs(os.path.join(td,'svsForAnnotation.csv'), os.path.join(td,'svsAnnotated.csv'), args.phylopPath, td)


# In[34]:

an = pd.read_csv(os.path.join(td,'svsAnnotated.csv'))


# In[35]:

# annotate SVs on each chromosome, using random forest trained on all other chroms, to avoid overfitting
an['path'] = 0
presentChroms = an['chrom'].value_counts().index.values
for chrm in presentChroms:
    rf = load('data/rfTrainedAllChromsExcept'+chrm+'.joblib')
    X = an[an['chrom'] == chrm][['DEL','numExonsFinal','phyloP', 'lowestExonRank', 'allSkippable','lowestExonsInGene', 'anyConstExon','pLIMax','loeufMin', 'cdsFracStartMin', 'cdsFracEndMax', 'cdsFracMax', 'pLI_max25_ID', 'loeuf_min25_ID','topExp','topUsage','maxStrength']].copy()
    an.loc[an['chrom'] == chrm,'path'] = rf.predict_proba(X)[:,1]


# In[36]:

an.set_index('OldID', inplace=True)


# Annotate vcf with StrVCTVRE pathogenicity scores

# In[42]:

if args.formatIn == 'vcf':
    print('\nwriting annotated VCF...\n')

    vcf = VCF(args.pathIn)
    vcf.add_info_to_header({'ID':'StrVCTVRE','Description':'pathogenicity score for structural variants','Type':'Float','Number':'1'})

    w = Writer(args.pathOut,vcf)

    count = 0
    for var in vcf:
        if var.INFO.get('END') and var.INFO.get('SVTYPE'):
            if count in an.index:
                var.INFO['StrVCTVRE'] = str(round(an.loc[count,'path'],3))
            elif count in validExon.index:
                var.INFO['StrVCTVRE'] = '1.0'
            elif count in validSVType.index:
                var.INFO['StrVCTVRE'] = 'not_exonic'
            elif count in validLength.index:
                var.INFO['StrVCTVRE'] = 'not_dup_or_del'    
            elif count in validStartEnd.index:
                var.INFO['StrVCTVRE'] = 'less_than_50bp'
            elif count in validChrom.index:
                var.INFO['StrVCTVRE'] = 'invalid_start_or_end'
            else:
                var.INFO['StrVCTVRE'] = 'not_valid_chrom'
            count += 1
        else:
            var.INFO['StrVCTVRE'] = 'missing_END_or_SVTYPE'
        w.write_record(var)
    w.close();
    vcf.close()
    
else:
    print('\nwriting annotated BED...\n')

    #bed = pd.read_csv(args.pathIn, sep='\t', names=['chrom','start','end','svtype'], header=None,dtype = {'chrom':str,'start':int,'end':int,'svtype':str})
    f = open(args.pathIn)
    outf = open(args.pathOut,'w')
    idx=0
    for row in [x.strip() for x in f.readlines()]:
        if idx in an.index:
            score = str(round(an.loc[idx,'path'],3))
        elif idx in validExon.index:
            score = '1.0'
        elif idx in validSVType.index:
            score = 'not_exonic'
        elif idx in validLength.index:
            score = 'not_dup_or_del'   
        elif idx in validStartEnd.index:
            score = 'less_than_50bp'
        elif idx in validChrom.index:
            score = 'invalid_start_or_end'
        else:
            score = 'not_valid_chrom'
        outf.write(row + '\t' + str(score) + '\n')
        idx += 1
    f.close()
    outf.close();


# delete temporary files

# In[43]:

shutil.rmtree(td)


# In[44]:

print('\nFinished\n')


# In[ ]:



