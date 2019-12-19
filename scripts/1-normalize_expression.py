"""Normalize Expression

Apply a quantile-normalization to the expression data, and standardize the
results so that each gene has a mean expression of zero.

Usage:
    normalize_expression.py [options] --rpkm=<file> --cnts=<file> --geno=<file> --attr=<file> --out=<dir>
    normalize_expression.py (-h | --help)

Options:
    --rpkm=<file>               File with RPKM.
    --cnts=<file>               File with reads count.
    --geno=<file>               VCF file with the genotypes.
    --attr=<file>               Sample attributes file.
    --out=<dir>                 Directory where the output will be saved.
    --rpkm_cutoff=<float>       Minimum RPKM. [default: 0.1]
    --cnts_cutoff=<int>         Minimum reads count. [default: 5]
    --sample_cutoff=<int>       Minimum number of samples with nonzero-expression. [default: 10]
    -h --help                   Print help

"""

# Adopted by Gao Wang from:
# https://github.com/broadinstitute/gtex-pipeline
# Originally authored by Francois Aguet
# Slightly modified by Federico Marotta (federico.marotta@edu.unito.it)

import numpy as np
import pandas as pd
import gzip
import subprocess
import scipy.stats as stats
import re, os
from docopt import docopt

def annotate_tissue_data(data, fsample):
    '''Save data to tissue specific tables'''
    samples = pd.read_csv(fsample, dtype=str, delimiter='\t', header=0)
    sample_dict = {}
    for row in samples[['SAMPID', 'SMTSD', 'SMAFRZE']].values:
        if row[2] == 'EXCLUDE':
            continue
        if row[1] not in sample_dict:
            sample_dict[row[1]] = []
        if row[0] in data.columns:
            sample_dict[row[1]].append(row[0])
    sample = dict((re.sub("[\W\d]+", "_", k.strip()).strip('_'), v) for k, v in sample_dict.items() if len(v))
    data = {k: data.loc[:, sample[k]] for k in sample}
    return data

def write_per_tissue_data(data, output_dir, suffix):
    # if os.path.isfile(output_dir):
    #    os.remove(output_dir)
    for k in data:
        # Strip specimen ID, keep only individual ID
        donor_ids = ['-'.join(i.split('-')[:2]) for i in data[k].columns]
        data[k].columns = donor_ids
        data[k].to_csv(output_dir + "/{k}.{suffix}".format(k = k, suffix = suffix), sep = '\t', compression = None)

def get_donors_from_vcf(vcfpath):
    """
    Extract donor IDs from VCF
    """
    with gzip.open(vcfpath) as vcf:
        for line in vcf:
            if line.decode()[:2]=='##': continue
            break
    return line.decode().strip().split('\t')[9:]

def read_gct(gct_file, donor_ids, dtype):
    """
    Load GCT as DataFrame
    First col of expression data is ENCODE gene name, 2nd col is HUGO name
    ======================================================================
    A more memory friendly version:

    head = pd.read_csv(fdata, skiprows = 2, sep = '\t', nrows = 1)
    dt = {'Description': str, 'Name': str}
    dt.update({x: dtype for x in head.columns if x not in dt})
    data = pd.read_csv(fdata, compression='gzip', skiprows=2,
                       index_col=0, header=0, dtype = dt, sep='\t').drop('Description', 1)

    """
    df = pd.read_csv(gct_file, sep='\t', skiprows=2, index_col=0)
    df.drop('Description', axis=1, inplace=True)
    df.index.name = 'GENE'
    return df[[i for i in df.columns if '-'.join(i.split('-')[:2]) in donor_ids]].astype(dtype, copy = True)

def normalize_quantiles(M, inplace=False):
    """
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")

    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003

    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    if not inplace:
        M = M.copy()

    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n

    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1

        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1

    if not inplace:
        return M

def inverse_quantile_normalization(M):
    """
    After quantile normalization of samples, standardize expression of each gene
    """
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q

def normalize_expression(expression_df, counts_df, expression_threshold=0.1, count_threshold=5, min_samples=10, dtype = np.float32):
    """
    Genes are thresholded based on the following expression rules:
      >=min_samples with >expression_threshold expression values
      >=min_samples with >count_threshold read counts
    """
    # donor_ids = ['-'.join(i.split('-')[:2]) for i in expression_df.columns]
    donor_ids = expression_df.columns

    # expression thresholds
    mask = ((np.sum(expression_df>expression_threshold,axis=1)>=min_samples) & (np.sum(counts_df>count_threshold,axis=1)>=min_samples)).values

    # apply normalization
    M = normalize_quantiles(expression_df.loc[mask].values, inplace=False)
    R = inverse_quantile_normalization(M)

    quant_std_df = pd.DataFrame(data=R, columns=donor_ids, index=expression_df.loc[mask].index, dtype = dtype)
    quant_df = pd.DataFrame(data=M, columns=donor_ids, index=expression_df.loc[mask].index, dtype = dtype)
    return quant_std_df, quant_df


# The real script starts here
args = docopt(__doc__, version='Normalize Expression 0.1.0')

print('Generating normalized expression files ... ', end='\n', flush=True)
donor_ids = get_donors_from_vcf(args["--geno"])
expression_df = read_gct(args["--rpkm"], donor_ids, np.float32)
counts_df = read_gct(args["--cnts"], donor_ids, np.uint32)
quant_std_df, quant_df = normalize_expression(expression_df, counts_df,
                                              expression_threshold=float(args["--rpkm_cutoff"]),
                                              count_threshold=int(args["--cnts_cutoff"]),
                                              min_samples=int(args["--sample_cutoff"]))
print('Save to TSV format, full matrix and per tissue data ...', end='\n', flush=True)
quant_std_per_tissue = annotate_tissue_data(quant_std_df, args["--attr"])
quant_per_tissue = annotate_tissue_data(quant_df, args["--attr"])
expression_per_tissue = annotate_tissue_data(expression_df, args["--attr"])
counts_per_tissue = annotate_tissue_data(counts_df, args["--attr"])
quant_df.to_csv(args["--out"] + "/normalized_expression.qnorm.tsv", sep = '\t', compression = None)
quant_std_df.to_csv(args["--out"] + "/normalized_expression.qnorm.std.tsv", sep = '\t', compression = None)
write_per_tissue_data(quant_per_tissue, args["--out"], "normalized_expression.qnorm.tsv")
write_per_tissue_data(quant_std_per_tissue, args["--out"], "normalized_expression.qnorm.std.tsv")
write_per_tissue_data(expression_per_tissue, args["--out"], "rpkm.tsv")
write_per_tissue_data(counts_per_tissue, args["--out"], "reads.tsv")
print('done.')
