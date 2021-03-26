import pandas as pd
import numpy as np
import asyncio
import pickle
import zipfile
import multiprocessing
import time
import argparse
import os
import warnings
from concurrent.futures import ProcessPoolExecutor


def merge_rs(chr_n, chunk, dat,  bpcol, trgt, reff):
    with zipfile.ZipFile(reff, 'r', compression=zipfile.ZIP_STORED) as zf:
        with zf.open('{}/{}'.format(chr_n, chunk), 'r') as f:
            rs = pickle.loads(f.read())
    dat[trgt] = dat[bpcol].map(rs)
    del rs
    nnasnps = np.where(~dat[trgt].isna())[0]
    trgtn = [n for n, i in enumerate(dat.columns) if i == trgt][0]
    dat.iloc[nnasnps, trgtn] = 'rs' + dat.iloc[nnasnps, trgtn].astype(np.int64).map(str)
    return dat


async def main(loop, sumstatsfile, outfile, chr_column, bp_column, rsid_column, reffile, threads):
    assert os.path.isfile(sumstatsfile), "Reference file not found"
    assert os.path.isfile(reffile), "Reference file not found"
    if os.path.isfile(outfile):
        warnings.warn("Output file already exists, this will be overwritten.")
    with zipfile.ZipFile(reffile, 'r', compression=zipfile.ZIP_STORED) as zf:
        with zf.open('_meta', 'r') as f:
            meta = pickle.loads(f.read())
    print("Reading input file   ", end='\r')
    dat = pd.read_csv(sumstatsfile, sep='\s+')
    assert chr_column in dat.columns, "Chromosome column not found"
    assert bp_column in dat.columns, "Basepair column not found"
    if rsid_column in dat.columns:
        warnings.warn("rsid column {} already exists, this will be replaced in the output.".format(rsid_column))
    chrin = dat[chr_column].astype(str) if 23 not in dat[chr_column].values else dat[chr_column].astype(str).replace({'23': 'X', '24': 'Y', '25': 'MT'})
    print("Splitting input      ", end='\r')
    userows = {i: dat.iloc[np.where(chrin.values == str(i))[0], :].copy() for i in meta['chromosomes']}
    del dat, chrin
    userows = {(i, j): dat.iloc[(dat[bp_column].values < ((j+1)*meta['chunksize'])) & (dat[bp_column].values >= (j*meta['chunksize'])), :].copy() for i, dat in userows.items() for j in range(meta['n_chunks'][i])}
    scriptargs = {k: [v, bp_column, rsid_column, reffile] for k, v in userows.items() if len(v) > 0}
    del userows
    print("Merging              ", end='\r')
    executor = ProcessPoolExecutor(max_workers=threads)
    data = await asyncio.gather(*(loop.run_in_executor(executor, merge_rs, *i, *v) for i, v in scriptargs.items()))
    print("Saving               ", end='\r')
    pd.concat(data, ignore_index=True, copy=False).to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Merge RSIDS to you Build38-based summary statistics")
    parser.add_argument("--sumstats", help="Input summary statistics, must be tab or space-delimited")
    parser.add_argument("--chr", help="Name of chromosome column in the input file", default='chr')
    parser.add_argument("--bp", help="Name of basepair position column in the input file", default='bp')
    parser.add_argument("--rsid", help="Optional name of creted/replaced RSID column", default='rsid')
    parser.add_argument("--out", help="Optional name of the output file (default = input+'_RS')", default=None)
    parser.add_argument("--ref", help="Path to .rsids reference file", default='rsids_GRCh38p7_build151.rsids')
    parser.add_argument("--threads", help="Number of threads to use (default is [detected threads]//2)", default=None)
    args = parser.parse_args()
    outf = args.sumstats+'_RS' if args.out is None else args.out

    n_threads = multiprocessing.cpu_count()-1 if args.threads is None else int(args.threads)
    if n_threads > (multiprocessing.cpu_count()-1):
        raise ValueError("Number of detected threads ({}) minus 1 is smaller than number of requested threads ({})".format(multiprocessing.cpu_count(), n_threads))
    start_t = time.time()
    lp = asyncio.get_event_loop()
    lp.run_until_complete(main(lp, args.sumstats, outf, args.chr, args.bp, args.rsid, args.ref, n_threads))
    dur = time.time() - start_t
    h, m = divmod(dur, 3600)
    m, s = divmod(m, 60)
    if h > 0:
        dur_fmt = '{:3d}:{:02d}:{:06.3f}'.format(int(h), int(m), s)
    else:
        dur_fmt = '{:02d}:{:06.3f}'.format(int(m), s) if m > 0 else '{:06.3f} seconds'.format(int(m), s)
    print('Script finished in '+dur_fmt)

