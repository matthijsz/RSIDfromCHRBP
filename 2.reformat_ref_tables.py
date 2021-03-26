import pandas as pd
import numpy as np
import pickle
import zipfile

chunksize = 1e6  # Size of individual pickle files in the .rsids file.
refname = 'rsids_GRCh38p7_build151'


def keep_last_build(x, axis):
    builds = x['dbSNPBuildID'].values
    idx = np.where(builds == builds.max())[0][0]
    return x.iloc[idx, 1:]


for i in list(range(1, 23)) + ['X', 'Y', 'MT']:
    print(i)
    dat = pd.read_csv("chr{}.INFO".format(i), sep='\t')
    buildat = pd.read_csv("buildchr{}.INFO".format(i), sep='\t')
    dat['dbSNPBuildID'] = buildat['dbSNPBuildID']
    del buildat
    dups = dat.duplicated(subset=['RS'], keep=False)
    dat_nodup = dat.loc[~dups, :].copy()
    dat_dup = dat.loc[dups, :].copy()
    if len(dat_dup) > 0:
        del dat
        print("Duplicates found in chr{}".format(i))
        # Remove SNPs with the same RSID in a single build
        dat_dup = dat_dup.loc[dat_dup.duplicated(subset=['RS', 'dbSNPBuildID'], keep=False), :].copy()
        dat_dup = dat_dup.groupby('RS').aggregate(axis=1, func=keep_last_build)
        dat = pd.concat([dat_nodup, dat_dup], ignore_index=True).reset_index()
        dat[['CHROM', 'POS', 'RS']].to_csv("chr{}.rsids".format(i), sep='\t', index=False)
    else:
        dat[['CHROM', 'POS', 'RS']].to_csv("chr{}.rsids".format(i), sep='\t', index=False)


rsids = {}
# Save .rsids Python file
nchunks = 0
chunk_layout = {'chunksize': chunksize, 'chromosomes': list(range(1, 23)) + ['X', 'Y', 'MT'],
                'n_chunks': {}}
with zipfile.ZipFile(refname + ".rsids", 'w', compression=zipfile.ZIP_STORED) as zf:
    for i in list(range(1, 23)) + ['X', 'Y', 'MT']:
        print('{}             '.format(i), end='\r')
        dat = pd.read_csv("chr{}.rsids".format(i), sep='\t')
        dat.index = dat['POS']
        rng = int((dat.index.max()//chunksize)+1)
        chunk_layout['n_chunks'][i] = rng
        for chunk in range(rng):
            with zf.open('{}/{}'.format(i, chunk), 'w') as f:
                print('{}:  {:4d}/{:4d}'.format(i, chunk, rng), end='\r')
                _ = f.write(pickle.dumps(dat.loc[(dat.index < (chunksize*(chunk+1))) & (dat.index >= chunksize*chunk), 'RS'].to_dict()))
            nchunks += 1
    with zf.open('_meta', 'w') as f:
        f.write(pickle.dumps(chunk_layout))

# Save .tsv.gz general file
for i in list(range(1, 23)) + ['X', 'Y', 'MT']:
    print(i)
    dat = pd.read_csv("chr{}.rsids".format(i), sep='\t')
    dat['RS'] = 'rs' + dat['RS'].map(str)
    if str(i) == "1":
        dat.to_csv(refname + ".tsv", sep='\t', index=False)
    else:
        dat.to_csv(refname + ".tsv", sep='\t', index=False, header=False, mode='a')