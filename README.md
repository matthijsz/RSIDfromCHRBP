## RSIDmerger

These scripts are to download and reformat reference VCF files to make it easy to merge RSID numbers to your GWAS results that only contain chromosome and basepair numbers.

Two files will be created:
 - a `.tsv.gz` file, this is a standard gzipped tab-seperated file with `CHR`, `BP` and `RS` columns.
 - a `.rsids` file, this is a special format designed to work with the accompanying `merge_rsids.py` for easy and (relatively) quick merging.

## Usage

####Generating reformatted reference files
 1. Run `1.download_ref_tables.sh` (bash script that uses vcftools) to download the VCF and extract the necessary columns. Change the link in that script, and the filename in the `vcftool` commands according to the reference builds you want.
 2. Run `2.reformat_ref_tables.py` to generate `.tsv.gz` and `.rsids` files. Change `refname` in this script depending on the builds downloaded
 
This takes some time but only needs to be done once per reference/build.
 
#### Merging RSIDs

##### Using merge_rsids.py

Run the script from command line as follows:
```
python merge_rs.py --sumstats [path to sumstatsfile] --chr [name of CHR column] --bp [name of BP column] --rsid [name of target RSID column] --out [name of output file] --ref [path to .rsids file] --threads [number of threads to use for multiprocessing]
```
The `--threads` argument will use all available threads (-1 for the main process) by default, and cannot be set to a number greater than this. Use this argument to decrease the number of threads used if you run into memory issues, or if you have a lot of threads (it doesn't scale well beyond 8 threads).
This will add a column to your file (name depends on input for `--rsid`) based on chromosome and basepair position.

Gzipped files can be used as both input and output, as `pandas` is used which automatically recognizes a `.gz` extension.

In my basic tests it took a little over 4 minutes to merge RSIDs from GRCh38-build151 (660M SNPs) to a file with about 3,5M SNPs.

##### Some other way

If you wish to merge RSIDs to your files in some other way I recommend you use the `.tsv.gz` files.
If you want to, you can also use the `.rsids` files as they are not some odd proprietary format (just a zipfile that holds pickle files, see the Python code).


 
