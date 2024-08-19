# Testing delimtools

### download data

```bash
# install
git clone https://github.com/boopsboops/ncbi-supermatrix.git
cd ncbi-supermatrix
Rscript -e "renv::restore()"
```

### run

```bash
# dry
scripts/download-sequences.R -c Chiloglanis -n 500 -x 2500 -b 1 -a false -d true
scripts/download-sequences.R -c Atopochilus -n 500 -x 2500 -b 1 -a false -d true

# real
scripts/download-sequences.R -c Chiloglanis -n 500 -x 2500 -b 500 -a false -d false
scripts/download-sequences.R -c Atopochilus -n 500 -x 2500 -b 10 -a true -d false

# cluster
scripts/clean-and-cluster.R -n 10 -c 0.6 -m 2

# pick
scripts/pick-clusters.R -c 3 -g cytb

# annotate
scripts/annotate-ncbi.R -t 1 -c fishbase

# filter
scripts/filter-species.R -n 3 -i true

# concat
scripts/align-trim-concatenate.R -p 0.2 -t 4 -i true

# tree
scripts/tree-search.R -m TN93+G -v false -e 10 -t 4

# plot
scripts/tree-plot.R -w 0.2 -h 0.8 -s 7

# clean
scripts/tidy-results-directory.R
```
