
# coding: utf-8

# # Assumes given a csv file with 'chrom', 'start', 'end' columns, convert from hg19 to hg38 or reverse

# In[1]:

import sys,os
import pandas as pd
import numpy as np


# In[2]:

def liftover(inp, work, hg19to38, outp, pathLiftover):
#     inp = sys.argv[1] # path to pandas dataframe csv that is currently in hg38
#     work = sys.argv[2] # path to working directory to save files
#     if sys.argv[3] == 'True':
#         hg19to38 = True
#     else:
#         hg19to38 = False
#     outp = sys.argv[4]

    # inp = '/data/andrewsharo-S/thesis/Aim1/AnnotSV/share/doc/AnnotSV/Example/test.bed.csv'
    # work = '/data/andrewsharo-S/thesis/Aim1'
    # hg19to38 = True
    # outp = '/data/andrewsharo-S/thesis/Aim1/AnnotSV/share/doc/AnnotSV/Example/test.bed.hg38.csv'

    original = work + '/toLO.bed'
    updated = work + '/fromLO.bed'
    error = work + '/fromLO.err'
    errorSansComments = work + '/fromLO.nocomments.err'
    errorLines = work+ '/toLO.bedtools.bed'
    originalFull = work + '/toMergeipynb.tsv'
    originalFullNoErrors = work + '/toMerge.bedtools.tsv'

    # Find chrom start and end columns, and make it so they are the first three columns in the front of the dataframe. This makes it easier to convert to a bed file

    df = pd.read_csv(inp)

    # first check that the chroms all have chr in front
    if sum(df['chrom'].astype(str).str.startswith('chr',na=False))/df.shape[0] < 0.5:
        df['chrom'] = 'chr' + df['chrom'].astype(str)

    inds = np.argwhere((df.columns.values=='chrom')|(df.columns.values=='start')|(df.columns.values=='end'))
    new_columns = np.delete(df.columns.values,inds)
    new_columns = np.append(['chrom', 'start', 'end'],new_columns)

    # For some unknown reason, liftover isn't able to map length 0 intervals. So add 1 to the end of any interval that is zero length

    df.loc[df['start'] == df['end'],'end'] = df['end'] + 1
    # ^ this looks sketchy, but works

    df[new_columns].to_csv(originalFull, index=False, sep='\t',header=False,na_rep='NaN')
    df_lo = df[['chrom', 'start','end']]
    df_lo.to_csv(original, index=False, sep='\t',header=False)

    # liftover to hg19 using local program

    if hg19to38:
        os.system(' '.join([pathLiftover, original, 'data/hg19ToHg38.over.chain.gz', updated, error]))
    else:
        os.system(' '.join([pathLiftover, original, 'data/hg38ToHg19.over.chain.gz', updated, error]))

    # Remove comments from the errors

    os.system(' '.join(['grep -v ^#', error, '>', errorSansComments]))

    # Remove errors from the original toLiftOver file, resulting in a file that will run cleanly through liftOver 

    os.system(' '.join(['bedtools intersect -v -f 1 -r -wa -a', original, '-b', errorSansComments, '>', errorLines]))

    # Remove errors from the original df, yielding a df that is error free

    os.system(' '.join(['bedtools intersect -v -f 1 -r -wa -a', originalFull, '-b', errorSansComments, '>', originalFullNoErrors]))

    # Read in error free original df, merge with converted df, save

    svOriginalNoErr = pd.read_csv(originalFullNoErrors, sep='\t',names=new_columns)
    if hg19to38:
        svUpdated = pd.read_csv(updated, sep='\t', names=['chrom_38','start_38','end_38'])
    else:
        svUpdated = pd.read_csv(updated, sep='\t', names=['chrom_37','start_37','end_37'])
    svUpdated.merge(svOriginalNoErr, left_index=True, right_index=True).to_csv(outp,index=False)


# In[ ]:




# In[ ]:



