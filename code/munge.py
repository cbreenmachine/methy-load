import os
import pandas as pd

# How to incorporate this script into GNU parallel;
# Need to pass the arguments through a bash script;
# Need to check what's been processed already ala extract script...
# Need to pass input file names (or a directory more likely) but let a bash script handle this

CHROMOSOMES = ['chr' + str(x) for x in range(1,23)]
CHROMOSOMES.append('chrX')
CHROMOSOMES.append('chrY')
ROOT = '../data/'
odir = '../data/cov_meth'
if not os.path.exists(odir):
    os.makedirs(odir)

def read_and_filter(file_path):
    sample = os.path.basename(file_path).replace('.tsv', '')
    # engine = 'python',
    df = pd.read_csv(file_path, sep="\t", usecols= ['chr', 'pos', 'methylated', 'unmethylated'], encoding="utf-8")

    df['coverage'] = df['methylated'] + df['unmethylated'] 
    df['sample'] = sample
    return df

def split_and_write(df):
    for c in CHROMOSOMES:
        outfile = os.path.join(odir, c + '.tsv')
        
        try:
            if not os.path.isfile(outfile):
                df[df['chr'] == c].to_csv(outfile, sep = '\t', header = 'column_names', index = False)
            else: # else it exists so append without writing the header
                #
                df[df['chr'] == c].to_csv(outfile, sep = '\t', mode = 'a', header = False, index = False)
        except UnicodeDecodeError:
            print(df['sample'][0] + " didn't write")


if __name__ == "__main__":

    all_files = []

    for root, subdir, files in os.walk(ROOT):
        if root.endswith("04-extract"):
            for ff in files:
                if ff.endswith(".tsv"):
                    all_files.append(os.path.join(root, ff))
    print(all_files)
    print("Working on " + str(len(all_files)) + " files")
    
    counter = 1
    for ff in all_files:
        print("Reading in " + ff)
        df = read_and_filter(ff)
        split_and_write(df)
        print("Finished {0} of {1}".format(counter, len(all_files)))
        counter += 1
