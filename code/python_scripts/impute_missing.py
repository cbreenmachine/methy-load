# want local mean with #of non-missing values in either direction
# local mean with max number of positions to search in either direction
# global mean
# global mean by category
# filter out cases where there is no methylation (0 or 1 methylated values) and lots of missing values
import pandas as pd
import argparse
import numpy as np
from sklearn.decomposition import IncrementalPCA
#

def read_data(file_path):
    df = pd.read_csv(file_path, sep = "\t")
    df_meth = (df.pivot_table(index=['pos'], columns=['sample'], values=['methylated']))
    df_cov = (df.pivot_table(index=['pos'], columns=['sample'], values=['coverage']))
    return df_meth, df_cov

def filter_too_many_nulls(df, num_na_cut=60, mean_cut=5):
    num_na = df.isnull().sum(axis=1)
    mean_val = df.sum(axis=1)

    ix = (num_na < num_na_cut) & (mean_val > mean_cut)
    return df[ix]


def impute_local_mean(df, group_ix, band_width=20):
    '''
    Args:
    df: a data frame with 
    group_ix = same length as number of columns in df

    '''
    return(None)
#https://scikit-learn.org/stable/auto_examples/decomposition/plot_incremental_pca.html


def run_pca(X, n_components=2, is_incremental=True):
    '''
    '''
    ipca = IncrementalPCA(n_components=n_components, batch_size=10000)
    X_ipca = ipca.fit_transform(X)
    return X_ipca


if __name__ == "main":
    # argparsing,...
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--ifile', default = '../../data/cov_meth/chr22_cov_meth.tsv', help='')
    parser.add_argument('--plot', action = 'store_true')
    args = parser.parse_args()

    df_meth, df_cov = read_data(args.ifile)
    
    #TODO: filter more sophisticated
    # Currently, filter columns with no nulls
    ix = df_meth.isnull().sum(axis=1) == 0
    df_meth_2 = df_meth[ix]
    
    X = df_meth_2.to_numpy().transpose()

    out = run_pca(X)
