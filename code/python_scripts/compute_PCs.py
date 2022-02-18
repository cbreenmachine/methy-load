#! /bin/python3
# compute_PCs.py takes
# RUN ON MASTODON--doesn't have memory issues there

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
    if 'methylation_estimate' not in df_raw:
        # point estimate for methylationg
        df_raw['methylation_estimate'] = df_raw['methylated'] / df_raw['coverage']
         
    df_raw.drop(columns=['methylated','chr', 'unmethylated'], inplace = True)
    df_raw = df_raw.astype({'sample': 'uint8', 'methylation_estimate': 'float32', 'coverage': 'uint8'})
   
    df = (df_raw.pivot_table(index=['pos'], columns=['sample'], values=['methylation_estimate', 'coverage']))
            
    return df

def filter_too_many_nulls(df):
    """Drops positions that are more than half nulls 
    
    """
    num_na = df.isnull().sum(axis=1)
    mean_val = df.sum(axis=1)

    ix = (num_na < num_na_cut) & (mean_val > mean_cut)
    return df[ix]


def impute_local_mean(df, group_ix, ws=50):
    # Designed with methylated reads and coverage in mind...
    
    '''imputes local mean by borrowin across groups
    
    Args:
    df: a data frame with 
    group_ix = same length as number of columns in df

    '''
    # Either mean then mean, or
    # Minimum periods := 
    mp = max(10, int(ws / 10))

    df.rolling(window = ws, min_periods = mp)


    return(None)

def run_pca(X, num_components = 2, is_incremental=True):
    '''computes principal components incremntally
    '''
    #TODO: allow for normal PCA
    ipca = IncrementalPCA(n_components = num_components, batch_size=10000)
    X_ipca = ipca.fit_transform(X)
    return X_ipca,  ipca.explained_variance_ratio_


if __name__ == "__main__":
    # argparsing,...
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--ifile', default = '../../data/cov-meth/chr22.tsv') #TODO: change to CSV in extract...
    parser.add_argument('--odir', default = '../../data/prin-comps-array-samples/')
    parser.add_argument('--filter_samples', action = 'store_true')
    parser.add_argument('--filter_file', default = '../../data/meta/array-samples.csv')
    args = parser.parse_args()
    
    
    if not os.path.exists(args.odir):
        os.makedirs(args.odir)

    if args.filter_samples:
        tmp = pd.read_csv(args.filter_file)['sample']
        filter_list = [(x) for x in tmp]

    df = read_data(args.ifile)
    X = df['methylation_estimate'].dropna(axis=0).transpose() # pos is index so gets dropped (no need to do anything else)    
    valid_samples = list(set(X.index).intersection(filter_list))
    X = X[X.index.isin(valid_samples)]

    # Don't drop nas for coverage--replace with zeros...
    Cov = df['coverage'].transpose()
    Cov = Cov[Cov.index.isin(valid_samples)]

    num_components = len(X) 
    print("Read in data frame, " + str(num_components) + " samples detected")
    
    # PCA step 
    pca_out, var_exp = run_pca(X, num_components = num_components)
    print("Computed PCs...")

    col_names = ['PC' + str(x) for x in range(1, num_components + 1)]
    pca_df = pd.DataFrame(pca_out, columns = col_names)
    
    # Add some other columns
    pca_df['sample'] = list(X.index)
    pca_df['num_null'] = list(Cov.isnull().sum(axis=1))
    pca_df['num_not_null'] = list(Cov.notnull().sum(axis=1))
    pca_df['mean_methylation'] = list(X.mean(axis=1, skipna=True))
    pca_df['mean_coverage'] = list(Cov.mean(axis=1, skipna=True))
    pca_df['median_coverage'] = list(Cov.median(axis=1, skipna=True))
    print("Computed summary statistics...")
    # file paths
    my_chr = os.path.basename(args.ifile).replace(".tsv", ".csv")
    ofile = os.path.join(args.odir, my_chr)
    pca_df.to_csv(ofile, index = False)
    print("Wrote out " + ofile)

    tmp = [var_exp, [i for i in range(1, num_components + 1)], [my_chr] * num_components]
    var_df = pd.DataFrame(tmp).T
    var_df.rename(columns = {0:'var_explained', 1:'PC', 2:'chrom'}, inplace = True)
    ofile = os.path.join(args.odir, 'var-explained-' + my_chr)

    #TODO: delete the if/else (deprecated) since we now write to individual (chr) files
    # if file does not exist write header 
    # https://stackoverflow.com/questions/30991541/pandas-write-csv-append-vs-write
    if not os.path.isfile(ofile):
        var_df.to_csv(ofile, header='column_names', index = False)
    else: # else it exists so append without writing the header
        var_df.to_csv(ofile, mode='a', header=False, index = False)



# Sandbox
Cov = df['coverage']

ix = [100, 101, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120, 122]

new_values = Cov[ix].mean(axis=1)

Cov[ix].fillna(new_values, axis=1)


def impute_with_group_mean(df, cols):
    new_values = df[cols].mean(axis = 1)

    for i, c in enumerate(df[cols]):
        df.iloc[:, i] = df.iloc[:, i].fillna(new_values)



Cov[cols] = Cov[cols].apply(lambda row: row.fillna(row.mean()), axis = 1)