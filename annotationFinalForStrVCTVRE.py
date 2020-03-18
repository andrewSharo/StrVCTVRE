
# coding: utf-8

# In[3]:

import pandas as pd
import pybedtools
import pyBigWig
import numpy as np
import os


# In[1]:

# Develop a function that takes the top N intervals from an SV, and averages those, instead of just doing max.

def topUsage(df,n):
    df['topUsage'] = df.sort_values('avgUsage',ascending=False).iloc[0:n,:]['avgUsage'].mean()
    df['topExp'] = df.sort_values('avgExp',ascending=False).iloc[0:n,:]['avgExp'].mean()
    return df


# In[2]:

# Define two functions that determine how far into the amino acids the SV begins

# assumes exon rank ordered.
def cdsRank(df):
    l = df.shape[0]
    dff = df.iloc[0,:] # get the first entry
    # first separate by strand
    if dff['strand'] == 1:
        if dff['start'] - dff['cStart'] < 0: # if sv begins before the first exon
            df['early'] = dff['cdsCount']
            return df
        else:
            df['early'] = dff['start'] - dff['cStart'] + dff['cdsCount']
            return df
    # reverse strand
    if dff['strand'] == -1:
        if dff['cStop'] - dff['stop'] < 0: # if sv begins before the first exon
            df['early'] = dff['cdsCount']
            return df
        else:
            df['early'] = dff['cStop'] - dff['stop'] + dff['cdsCount']
            return df
        
# assumes reverse exon rank ordered.
def cdsEnd(df):
    l = df.shape[0]
    dff = df.iloc[0,:]
    # first separate by strand
    if dff['strand'] == 1:
        if dff['stop'] > dff['cStop']: # if sv begins before the first exon
            df['late'] = dff['cdsCount'] + dff['size']
            return df
        else:
            df['late'] = dff['stop'] - dff['cStart'] + dff['cdsCount']
            return df
    # reverse strand
    if dff['strand'] == -1:
        if dff['start'] < dff['cStart']: # if sv ends after the last overlapped exon
            df['late'] = dff['cdsCount'] + dff['size']
            return df
        else:
            df['late'] = dff['cStop'] - dff['start'] + dff['cdsCount']
            return df


# In[2]:

# for def:
def annotateSVs(inpath, outpath, phylopPath, tempdir):

    # read csv file into dataframe
    
    df = pd.read_csv(inpath)
    
    # Do all exon-level and gene-level features

    exons = pybedtools.BedTool('data/exons_Appris_featurized_transcript_Chr1-Y_loeuf.sorted.bed')
    df['ID'] = 'sv' + pd.Series(df.index.values).apply(str)
    df[['chrom','start','end','ID']].to_csv(os.path.join(tempdir,'df.bed'),sep='\t', index=False,header=False)
    a = pybedtools.BedTool(os.path.join(tempdir,'df.bed'))
    b = a.intersect(exons, wa=True, wb=True).saveas(os.path.join(tempdir,'dfExonOverlap.bed'))
    exonOverlap = pd.read_csv(os.path.join(tempdir,'dfExonOverlap.bed'), sep='\t', header=None, 
                              names=['chrom', 'start', 'stop', 'ID', 'eChrom', 'eStart', 'eStop', 'gene', 'exonRank', 'skippable', 'exonsInGene', 'const','pLI','loeuf'])

    exonOverlap['numExonsFinal'] = exonOverlap.groupby(by='ID').eStart.transform('size')
    exonOverlap['allSkippable'] = exonOverlap.groupby(by='ID').skippable.transform(lambda x: all(x))
    exonOverlap['lowestExonRank'] = exonOverlap.groupby(by='ID').exonRank.transform('min')
    exonOverlap['lowestExonsInGene'] = exonOverlap.groupby(by='ID').exonsInGene.transform('min')
    exonOverlap['anyConstExon'] = exonOverlap.groupby(by='ID').const.transform('max')
    exonOverlap['pLIMax'] = exonOverlap.groupby(by='ID').pLI.transform('max')
    exonOverlap['loeufMin'] = exonOverlap.groupby(by='ID').loeuf.transform('min')

    exonOverlap.drop_duplicates(subset='ID', inplace=True)
    df = df.merge(exonOverlap[['ID', 'numExonsFinal', 'allSkippable', 'lowestExonRank', 'lowestExonsInGene', 'anyConstExon','pLIMax','loeufMin']], how='left', on='ID')

    # numExons = the total number of exons that an SV overlaps, across all genes

    # allSkippable = 1 if all exons overlapped start and end in same phase, 0 otherwise

    # lowestExonRank = the minimum rank of all exons overlapped

    # lowestExonsInGene = the number of exons in the gene overlapped, minimum if multiple genes

    # anyConstExon = 1 if any exon overlapped is constitutive, 0 otherwise

    # pLIMax = the maximum pLI of all overlapped genes (high pLI is more intolerant)

    # loeufMin = the minimum LOEUF of all overlapped genes (low loeuf is more intolerant)

    # Calculate conservation feature using phyloP as average of top 400 most conserved position 

    size = 400

    with open("data/hg38chromsizes.tsv") as f:
        chrms = dict((k, v) for k,v in (line.split() for line in f))

    consBW = pyBigWig.open(phylopPath)
    # get phyloP value for each position in the SV
    x = []
    for i in range(df.shape[0]):
        if int(df.loc[i,'end']) - int(df.loc[i,'start']) > 1000000:
            x.append(np.array([15.0]))
        else:
            try:
                x.append(np.nan_to_num(np.array(consBW.values(df.loc[i,'chrom'], int(df.loc[i,'start']), int(df.loc[i,'end'])))))
            except:
                print(df.loc[i,'start'])
                x.append(np.array([0.5]))
    x = np.asarray(x)
    # get the mean of the top 100 most conserved positions
    cons = [np.mean(y[np.argsort(y)[-size:]]) for y in x]
    df['phyloP'] = pd.Series(cons)

    # Add TAD features

    tads = pybedtools.BedTool('data/rep12tadsMergedhg38.bed')
    df[['chrom','start','end','ID']].to_csv(os.path.join(tempdir,'df.bed'),sep='\t', index=False,header=False)
    a = pybedtools.BedTool(os.path.join(tempdir,'df.bed'))
    b = a.intersect(tads, wa=True, wb=True).saveas(os.path.join(tempdir,'dfTadOverlap.bed'))
    tadOverlap = pd.read_csv(os.path.join(tempdir,'dfTadOverlap.bed'), sep='\t', header=None, 
                              names=['chrom', 'start', 'stop', 'ID', 'tChrom', 'tStart', 'tStop', 'strength'])

    tadOverlap['maxStrength'] = tadOverlap.groupby(by='ID').strength.transform('max')
    tadOverlap.drop_duplicates(subset='ID', inplace=True)
    df = df.merge(tadOverlap[['ID', 'maxStrength']], how='left', on='ID')
    df['maxStrength'].fillna(value=0, inplace=True)

    ## Add amino acid features

    cds = pybedtools.BedTool('data/exons_CDS_Chr1-Y.sorted.bed')
    df[['chrom','start','end','ID']].to_csv(os.path.join(tempdir,'df.bed'),sep='\t', index=False,header=False)
    a = pybedtools.BedTool(os.path.join(tempdir,'df.bed'))
    b = a.intersect(cds, wa=True, wb=True).saveas(os.path.join(tempdir,'dfCDSOverlap.bed'))
    cdsOverlap = pd.read_csv(os.path.join(tempdir,'dfCDSOverlap.bed'), sep='\t', header=None, 
                              names=['chrom', 'start', 'stop', 'ID', 'cChrom', 'cStart', 'cStop', 'CDSLength', 'size', 'exonRank', 'strand','gene', 'cdsCount', 'pLI','loeuf'])
    
    # use if statement to address possible scenario in which all given variants are in UTR and don't overlap a CDS
    if cdsOverlap.shape[0] != 0:
        # apply above functions to the SVs that were previously intersected with coding exons
        out = cdsOverlap.sort_values('exonRank').groupby(['ID','gene']).apply(cdsRank)
        out = out.sort_values('exonRank', ascending=False).groupby(['ID','gene']).apply(cdsEnd)

        # get shape, but don't drop duplicates (not in place)
        out.drop_duplicates(subset='ID').shape[0]

        # Featurize above information into features normalized by cds length

        out['cdsFracStart'] = out['early']/out['CDSLength']
        out['cdsFracEnd'] = out['late']/out['CDSLength']
        out['cdsFrac'] = (out['late'] - out['early'])/out['CDSLength'] 

        # This is an experimental feature, which gives the max pLI and loeuf of the genes which are signficantly disrupted by the SV.

        out['pLI_max25'] = out[(out['cdsFracStart'] == 0) | (out['cdsFrac'] > 0.25)].groupby('ID')['pLI'].transform('max')
        #out['pLI_max25'].fillna(value=0, inplace=True)
        out['loeuf_min25']= out[(out['cdsFracStart'] == 0) | (out['cdsFrac'] > 0.25)].groupby('ID')['loeuf'].transform('min')
        #out['loeuf_min25'].fillna(value=0, inplace=True)

        # but we now need to fill in all the cells with the max loeuf_max25 in their ID
        out['pLI_max25_ID'] = out.groupby('ID')['pLI_max25'].transform('max')
        out['loeuf_min25_ID'] = out.groupby('ID')['loeuf_min25'].transform('max')

        out['cdsFracMax'] = out.groupby('ID')['cdsFrac'].transform('max')
        out['cdsFracStartMin'] = out.groupby('ID')['cdsFracStart'].transform('min')
        out['cdsFracEndMax'] = out.groupby('ID')['cdsFracEnd'].transform('max')

        out.drop_duplicates(subset='ID', inplace=True)
        
        final = df.merge(out[['ID', 'cdsFracStartMin', 'cdsFracEndMax', 'cdsFracMax', 'pLI_max25_ID', 'loeuf_min25_ID']], how='left')
    else:
        final = df.copy()
        
        final['cdsFracStartMin'] = float('NaN')
        final['cdsFracEndMax'] = float('NaN')
        final['cdsFracMax'] = float('NaN')
        final['pLI_max25_ID'] = float('NaN')
        final['loeuf_min25_ID'] = float('NaN')


    final['cdsFracStartMin'].fillna(value=2, inplace=True)
    final['cdsFracEndMax'].fillna(value=-1, inplace=True)
    final['cdsFracMax'].fillna(value=-1, inplace=True)
    final['pLI_max25_ID'].fillna(value=-1, inplace=True)
    final['loeuf_min25_ID'].fillna(value=3, inplace=True)

    # Add exon inclusion features

    usage = pybedtools.BedTool('data/summary_exon_usage_hg38.sorted.bed')
    final[['chrom','start','end','ID']].to_csv(os.path.join(tempdir,'df.bed'),sep='\t', index=False,header=False)
    a = pybedtools.BedTool(os.path.join(tempdir,'df.bed'))
    b = a.intersect(usage, wa=True, wb=True).saveas(os.path.join(tempdir,'dfUsageOverlap.bed'))
    usageOverlap = pd.read_csv(os.path.join(tempdir,'dfUsageOverlap.bed'), sep='\t', header=None, 
                              names=['chrom', 'start', 'stop', 'ID', 'uChrom', 'uStart', 'uStop', 'avgUsage', 'avgExp'])
    out = usageOverlap.groupby('ID').apply(topUsage,n=size)
    out.drop_duplicates(subset='ID', inplace=True)

    final2 = final.merge(out[['ID', 'topUsage', 'topExp']], how='left', on='ID')
    final2['topExp'].fillna(value=final2['topExp'].median(), inplace=True)
    final2['topUsage'].fillna(value=final2['topUsage'].median(), inplace=True)

    final2.to_csv(outpath,index=False)

