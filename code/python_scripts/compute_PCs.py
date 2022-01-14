#! /bin/python3
# impute_missing.py
# AS it stands, this is more of a PCA script than anything else...


# want local mean with #of non-missing values in either direction
# local mean with max number of positions to search in either direction
# global mean
# global mean by category
# filter out cases where there is no methylation (0 or 1 methylated values) and lots of missing values
import pandas as pd
import argparse
import os
from sklearn.decomposition import IncrementalPCA


def read_data(file_path):
    '''
    Parameters
    ----------
        file_path : str
            relative path to one chromosome
    
    Returns
    -------

    '''

    df_raw = pd.read_csv(file_path, sep = "\t")
    df_raw['methylation_estimate'] = df_raw['methylated'] / df_raw['coverage'] # point estimate for methylationg
    df_raw.drop(columns=['methylated','chr', 'unmethylated'], inplace = True)
    df_raw = df_raw.astype({'sample': 'uint8', 'methylation_estimate': 'float32', 'coverage': 'uint8'})
    # TODO: Delet these lines if this works...
    #df_meth = (df.pivot_table(index=['pos'], columns=['sample'], values=['methylated']))
    #df_cov = (df.pivot_table(index=['pos'], columns=['sample'], values=['coverage']))
    df = (df_raw.pivot_table(index=['pos'], columns=['sample'], values=['methylation_estimate', 'coverage']))
    return df

def filter_too_many_nulls(df, num_na_cut=60, mean_cut=5):
    """Drops positions that are (almost) all nulls 
    
    """
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


def run_pca(X, num_components = 2, is_incremental=True):
    '''
    '''
    #TODO: allow for normal PCA
    ipca = IncrementalPCA(n_components = num_components, batch_size=10000)
    X_ipca = ipca.fit_transform(X)
    return X_ipca,  ipca.explained_variance_ratio_


if __name__ == "__main__":
    # argparsing,...
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--ifile', default = '../../data/cov_meth/chr22_cov_meth.tsv')
    parser.add_argument('--odir', default = '../../data/prin_comps/')
    parser.add_argument('--plot', action = 'store_true')
    parser.add_argument('--num_components', default = 10)
    args = parser.parse_args()
    
    if not os.path.exists(args.odir):
        os.makedirs(args.odir)

    df = read_data(args.ifile)
    
    #TODO: filter more sophisticated
    # Currently, filter columns with no nulls
    ix = df['methylation_estimate'].isnull().sum(axis=1) == 0
    df_2 = df['methylation_estimate'][ix] # subset
    X = df_2.to_numpy().transpose() # pos is index so gets dropped (no need to do anything else)    

    # PCA step 
    # TODO: Run normal PCA or incremental??
    pca_out, var_exp = run_pca(X, num_components = args.num_components)
    # 



    col_names = ['PC' + str(x) for x in range(1, args.num_components + 1)]
    pca_df = pd.DataFrame(pca_out, columns = col_names)
    
    # Add some other columns
    pca_df['sample'] = list(df.columns.levels[1])
    pca_df['num_null'] = list(df['methylation_estimate'].isnull().sum(axis=0))
    pca_df['num_not_null'] = list(df['methylation_estimate'].notnull().sum(axis=0))
    pca_df['mean_methylation'] = list(df['methylation_estimate'].mean(axis=0, skipna=True))
    pca_df['mean_coverage'] = list(df['coverage'].mean(axis=0, skipna=True))
    pca_df['median_coverage'] = list(df['coverage'].median(axis=0, skipna=True))

    # file paths
    my_chr = os.path.basename(args.ifile).replace("_cov_meth.tsv", "")
    ofile = os.path.join(args.odir, 'PC_' + my_chr + '.tsv')
    pca_df.to_csv(ofile, sep = '\t', index = False)
    print("Wrote out " + ofile)

    tmp = [var_exp, [i for i in range(1, args.num_components + 1)], [my_chr] * args.num_components]
    var_df = pd.DataFrame(tmp).T
    var_df.rename(columns = {0:'var_explained', 1:'PC', 2:'chrom'}, inplace = True)
    #var_df
    ofile = os.path.join(args.odir, 'var_explained.tsv')

    # if file does not exist write header 
    # https://stackoverflow.com/questions/30991541/pandas-write-csv-append-vs-write
    if not os.path.isfile(ofile):
        var_df.to_csv(ofile, header='column_names', sep = '\t', index = False)
    else: # else it exists so append without writing the header
        var_df.to_csv(ofile, mode='a', header=False, sep = '\t', index = False)

