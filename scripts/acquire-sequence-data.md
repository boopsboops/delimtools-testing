# How to access testing data for use with delimtools


### Run ncbi-supermatrix to obtain a COI dataset of Geophagus from Ximenes et al. (2021)
##### !!! RUN IN BASH !!!

```bash
# need to be in 'delimtools-testing' 
cd ncbi-supermatrix
scripts/download-sequences.R -c Geophagus -n 500 -x 2500 -b 1 -a false -d true
scripts/download-sequences.R -c Geophagus -n 500 -x 2500 -b 500 -a false -d false
scripts/clean-and-cluster.R -n 10 -c 0.6 -m 2
scripts/pick-clusters.R -c 10 -g coi
scripts/annotate-ncbi.R -t 3 -c fishbase
scripts/filter-species.R -n 3 -i true
scripts/align-trim-concatenate.R -p 0.2 -t 4 -i true
scripts/tree-search.R -m TN93+G -v false -e 10 -t 4
scripts/tree-plot.R -w 0.2 -h 0.8 -s 7
scripts/tidy-results-directory.R
```


### Subset samples for ingroup (Geophagus sensu stricto) and generate input data
##### !!! RUN IN R !!! ###

```r
# install renv for testing
#renv::install(here::here("../delimtools"))
#renv::install("legalLab/delimtools")
#renv::install(here::here(getwd(),'software/bGMYC_1.0.3.tar.gz'))
source(here::here("scripts/load-libs.R"))

# create a temporary directory
today.dir <- glue::glue('Results_{Sys.Date()}')
today.path <- here::here("temp",today.dir)
if(!dir.exists(today.path)) {dir.create(today.path,recursive=TRUE)}

# copy files from ncbi-supermatrix
file.copy(here::here("ncbi-supermatrix/temp/Results_2024-08-27/alignments/coi.aligned.trimmed.fasta"), here::here(today.path,"coi.geophagus.fasta"))
file.copy(here::here("ncbi-supermatrix/temp/Results_2024-08-27/trees/coi.aligned.trimmed.fasta.raxml.bestTree"), here::here(today.path,"coi.geophagus.nwk"))
file.copy(here::here("ncbi-supermatrix/temp/Results_2024-08-27/metadata/ncbi-clean.csv"), here::here(today.path,"coi.geophagus.csv"))

# read fasta
coi.geophagus.fa <- ape::read.FASTA(here::here(today.path,"coi.geophagus.fasta"))
# read tree
coi.geophagus.raxml.tr <- ape::read.tree(here::here(today.path,"coi.geophagus.nwk"))
# read data 
coi.geophagus.df <- readr::read_csv(here::here(today.path,"coi.geophagus.csv"),show_col_types=FALSE)

# get ingroup and outgroups
ingroup.tips <- coi.geophagus.df |> 
    dplyr::filter(grepl("harreri|winemilleri",scientificName)) |> 
    dplyr::slice_head(n=1,by=scientificName) |> 
    dplyr::pull(gbAccession)

# get one outgroup (any G. brasiliensis)
outgroup <- coi.geophagus.df |> 
    dplyr::filter(grepl("brasiliensis",scientificName)) |> 
    dplyr::slice_head(n=1) |> 
    dplyr::pull(gbAccession)

# get all ingroup tips
coi.geophagus.raxml.tr.rooted <- coi.geophagus.raxml.tr |> ape::root(outgroup=outgroup,resolve.root=TRUE)
ingroup.node <- coi.geophagus.raxml.tr.rooted |> ape::getMRCA(ingroup.tips)
ingroup <- ape::extract.clade(coi.geophagus.raxml.tr.rooted,node=ingroup.node)$tip.label

# join ingroup and outgroup
geo.gb <- c(ingroup,outgroup)

# subset ingroup from fasta
coi.geophagus.ingroup.fa <- coi.geophagus.fa[geo.gb]

# collapse haplotypes
delimtools::haplotype_tbl(coi.geophagus.ingroup.fa,verbose=FALSE) |> print(n=Inf)
haps.collapsed <- delimtools::hap_collapse(coi.geophagus.ingroup.fa,collapseSubstrings=TRUE,clean=TRUE) |> names()

# subset and write out
coi.geophagus.haps.fa <- coi.geophagus.ingroup.fa[haps.collapsed]
coi.geophagus.haps.fa |> ape::write.FASTA(here(today.path,"coi.geophagus.haps.fasta"))

# write out df
coi.geophagus.df |> 
    filter(gbAccession %in% haps.collapsed) |> 
    readr::write_csv(here::here(today.path,"coi.geophagus.haps.csv"))
```


### MAKE PHYLOGENETIC TREES USING BEAST AND RAXML
##### !!! RUN IN BASH !!!

```bash
# need to be in 'delimtools-testing'
# find most recent tempdir
tmpdir=$(ls -d temp/*/ | sort -r | head -n 1)
echo $tmpdir
cd $tmpdir
# beauti coi.geophagus.haps.fasta
# set up beauti for two .xml files using strict clock, constant coalescent, TN93+G model, chain length 10000000, log every 16000
cp ../..assets/coi.geophagus.haps.run1.xml coi.geophagus.haps.run1.xml
cp ../..assets/coi.geophagus.haps.run2.xml coi.geophagus.haps.run2.xml 

# run beast
beast -beagle_auto -overwrite -seed 42 coi.geophagus.haps.run1.xml 
beast -beagle_auto -overwrite -seed 24 coi.geophagus.haps.run2.xml 
tracer coi.geophagus.haps.run1.log coi.geophagus.haps.run2.log
logcombiner -trees -burnin 2016000 coi.geophagus.haps.run1.trees coi.geophagus.haps.run2.trees coi.geophagus.haps.run1+2.trees
grep -c "tree STATE_" coi.geophagus.haps.run1+2.trees
treeannotator -burninTrees 0 -heights ca coi.geophagus.haps.run1+2.trees coi.geophagus.haps.beast.tre
# cp coi.geophagus.haps.run1+2.trees ../..assets/coi.geophagus.haps.run1+2.trees

# run raxml on beast tree
# run in R first to convert tree to newick
#ape::read.nexus(here(today.path,"coi.geophagus.haps.beast.tre")) |> write.tree(here(today.path,"coi.geophagus.haps.beast.tre.nwk"))
raxml-ng --evaluate --threads auto{} --tree coi.geophagus.haps.beast.tre.nwk --lh-epsilon 0.1 --redo --seed 42 --outgroup MH538063.1 --model TN93+G --msa coi.geophagus.haps.fasta
```


### PROCESS PHYLOGENETIC TREES, DROP OUTGROUP AND SAVE
##### !!! RUN IN R !!!

```r
# get temp dir
today.dir <- glue::glue('Results_{Sys.Date()}')
today.path <- here::here("temp",today.dir)
if(!dir.exists(today.path)) {dir.create(today.path,recursive=TRUE)}

# read tree and drop outgroup
coi.geophagus.haps.raxml.tr <- ape::read.tree(here::here(today.path,"coi.geophagus.haps.fasta.raxml.bestTree")) |> ape::drop.tip("MH538063.1")

# read beast tree and drop outgroup
coi.geophagus.haps.beast.tr <- treeio::read.beast(here::here(today.path,"coi.geophagus.haps.beast.tre")) |> treeio::drop.tip("MH538063.1")

# read data 
coi.geophagus.haps.df <- readr::read_csv(here::here(today.path,"coi.geophagus.haps.csv"),show_col_types=FALSE) |> filter(gbAccession!="MH538063.1")

# read dna 
coi.geophagus.haps.fa <- ape::read.FASTA(here::here(today.path,"coi.geophagus.haps.fasta"))
coi.geophagus.haps.fa <- coi.geophagus.haps.fa[which(names(coi.geophagus.haps.fa)!="MH538063.1")]

# write out
coi.geophagus.haps.raxml.tr |> ape::write.tree(here("assets/coi.geophagus.haps.raxml.nwk"))
coi.geophagus.haps.beast.tr |> treeio::write.beast(here("assets/coi.geophagus.haps.beast.tre"))
coi.geophagus.haps.df |> readr::write_csv(here("assets/coi.geophagus.haps.csv"))
coi.geophagus.haps.fa |> ape::write.FASTA(here("assets/coi.geophagus.haps.fasta"))
```
