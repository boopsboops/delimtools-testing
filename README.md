# Testing delimtools

### set up
```bash
git clone https://github.com/boopsboops/delimtools-testing.git
cd delimtools-testing
mkdir temp software
cp assets/splits_1.0-20.tar.gz software/splits_1.0-20.tar.gz
cp assets/bGMYC_1.0.2.tar.gz software/bGMYC_1.0.2.tar.gz
```

### install software

```bash
# install raxml-ng
cd software
wget https://github.com/amkozlov/raxml-ng/releases/download/1.2.2/raxml-ng_v1.2.2_linux_x86_64.zip
unzip raxml-ng_v1.2.2_linux_x86_64.zip -d raxml-ng
echo "export PATH=$(pwd)/raxml-ng:\$PATH" >> ~/.bashrc
exec "$SHELL"
```

```bash
# install mptp
sudo apt-get install libgsl0-dev flex bison autotools-dev autoconf
git clone https://github.com/Pas-Kapli/mptp.git
cd mptp
git checkout v0.2.5
./autogen.sh
./configure
make
echo "export PATH=$(pwd)/bin:\$PATH" >> ~/.bashrc
exec "$SHELL"
cd ..
```

```bash
wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0-beta4/BEAST_X_v10.5.0-beta4.tgz
tar -xzvf BEAST_X_v10.5.0-beta4.tgz
echo "export PATH=$(pwd)/BEASTv10.5.0/bin:\$PATH" >> ~/.bashrc
exec "$SHELL"
cd ..
```

```bash
wget https://bioinfo.mnhn.fr/abi/public/asap/last.tgz
tar -xzvf last.tgz
cd ASAP
make
mkdir bin
mv asap bin/asap
echo "export PATH=$(pwd)/bin:\$PATH" >> ~/.bashrc
exec "$SHELL"
cd ../..
```


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
scripts/download-sequences.R -c Euchilichthys -n 500 -x 2500 -b 1 -a false -d true
# real
scripts/download-sequences.R -c Chiloglanis -n 500 -x 2500 -b 500 -a false -d false
scripts/download-sequences.R -c Atopochilus -n 500 -x 2500 -b 10 -a true -d false
scripts/download-sequences.R -c Euchilichthys -n 500 -x 2500 -b 20 -a true -d false
# cluster
scripts/clean-and-cluster.R -n 10 -c 0.6 -m 2
# pick
scripts/pick-clusters.R -c 3 -g cytb
# annotate
scripts/annotate-ncbi.R -t 2 -c fishbase
# filter
scripts/filter-species.R -n 3 -i true
# concat
scripts/align-trim-concatenate.R -p 0.2 -t 4 -i true
# quick tree
scripts/tree-search.R -m TN93+G -v false -e 10 -t 4
# plot
scripts/tree-plot.R -w 0.2 -h 0.8 -s 7
# clean
scripts/tidy-results-directory.R

# copy
cp temp/Results_2024-08-20/alignments/cytb.aligned.trimmed.fasta ../assets/cytb.chiloglanis.fasta
cp temp/Results_2024-08-20/trees/cytb.aligned.trimmed.fasta.raxml.bestTree ../assets/cytb.chiloglanis.nwk
cp temp/Results_2024-08-20/metadata/ncbi-clean.csv ../assets/cytb.chiloglanis.csv
```

### install R packages for delim tools

```bash
Rscript -e "renv::restore()"
# need to sort out bGMYC and GMYC package installs
# wget https://nreid.github.io/assets/bGMYC_1.0.2.tar.gz
# wget http://download.r-forge.r-project.org/src/contrib/splits_1.0-20.tar.gz
Rscript -e "renv::install(here::here(getwd(),'software/splits_1.0-20.tar.gz'))"
Rscript -e "renv::install(here::here(getwd(),'software/bGMYC_1.0.3.tar.gz'))"
```



### generate trees

```r
# install renv for testing
#renv::install(here::here("../delimtools"))
#renv::install("legalLab/delimtools")
#renv::install(here::here(getwd(),'software/bGMYC_1.0.3.tar.gz'))
library("here")
library("glue")
library("tidyverse")
library("ape")
library("spider")
library("splits")
library("bGMYC")
library("delimtools")

# 
rm(list=ls())

# read tree
cytb.chiloglanis.raxml.tr <- ape::read.tree(here("assets/cytb.chiloglanis.nwk"))
cytb.chiloglanis.beast.tr <- ape::read.nexus(here("temp/beast/cytb.chiloglanis.beast.tre"))

# read data 
cytb.chiloglanis.df <- read_csv(here("assets/cytb.chiloglanis.csv"),show_col_types=FALSE)
#cytb.chiloglanis.df |> glimpse()
# read fasta
cytb.chiloglanis.fa <- ape::read.FASTA(here("assets/cytb.chiloglanis.fasta"))

# collapse haplotypes
# hap table
delimtools::haplotype_tbl(cytb.chiloglanis.fa,verbose=FALSE) |> print(n=Inf)
haps.collapsed <- delimtools::hap_collapse(cytb.chiloglanis.fa,collapseSubstrings=TRUE,clean=TRUE) |> names()

# root tree
outgroups <- cytb.chiloglanis.df |> 
    filter(genus!="Chiloglanis") |> 
    pull(gbAccession)

# undescribed and >5
undescribed <- cytb.chiloglanis.df |> 
    filter(grepl("sp\\.|cf\\.|aff\\.",scientificName)) |> 
    pull(gbAccession)

# get root
root.node <- ape::getMRCA(cytb.chiloglanis.raxml.tr,tip=outgroups)

# remove > 5
set.seed(42)
sliced <- cytb.chiloglanis.df |>
    group_by(scientificName) |>
    slice_sample(n=5) |>
    ungroup() |> 
    pull(gbAccession)
excess <- cytb.chiloglanis.df |> filter(!gbAccession %in% sliced) |> pull(gbAccession)

# join drop
#drops <- unique(c(setdiff(names(cytb.chiloglanis.fa),haps.collapsed),outgroups))
drops <- unique(c(setdiff(names(cytb.chiloglanis.fa),haps.collapsed),outgroups,undescribed,excess))

# root and drop outgroups and duplicated haps
cytb.chiloglanis.raxml.tr.root <- cytb.chiloglanis.raxml.tr |> 
    ape::root(node=root.node,resolve.root=TRUE) |> 
    ape::drop.tip(drops)

# root and drop outgroups and duplicated haps
cytb.chiloglanis.beast.tr.root <- cytb.chiloglanis.beast.tr |> 
    ape::drop.tip(drops)

# save out
# make date dir for results
today.dir <- glue('Results_{Sys.Date()}')
today.path <- here("temp",today.dir)
if(!dir.exists(today.path)) {dir.create(today.path,recursive=TRUE)}

# drop from fasta and write
cytb.chiloglanis.fa.haps <- cytb.chiloglanis.fa[which(!names(cytb.chiloglanis.fa) %in% drops)]

# save
cytb.chiloglanis.fa.haps |> ape::write.FASTA(here(today.path,"cytb.chiloglanis.haps.fasta"))
raxml.tr.path <- here(today.path,"cytb.chiloglanis.rooted.nwk")
beast.tr.path <- here(today.path,"cytb.chiloglanis.rooted.tre")
cytb.chiloglanis.raxml.tr.root |> ape::write.tree(file=raxml.tr.path)
cytb.chiloglanis.beast.tr.root |> ape::write.nexus(file=beast.tr.path)

#cytb.chiloglanis.df |> filter(gbAccession %in% names(cytb.chiloglanis.fa.haps))
#cytb.chiloglanis.df |> filter(gbAccession %in% cytb.chiloglanis.raxml.tr.root$tip.label)
#cytb.chiloglanis.df |> filter(gbAccession %in% cytb.chiloglanis.beast.tr.root$tip.label)

```

```bash
beauti assets/cytb.chiloglanis.fasta
mkdir temp/beast
cp assets/cytb.chiloglanis.run.1.xml temp/beast/cytb.chiloglanis.run.1.xml 
cp assets/cytb.chiloglanis.run.2.xml temp/beast/cytb.chiloglanis.run.2.xml 

# run beast
cd temp/beast
beast -beagle_auto -overwrite -seed 42 cytb.chiloglanis.run.1.xml 
beast -beagle_auto -overwrite -seed 41 cytb.chiloglanis.run.2.xml 
tracer cytb.chiloglanis.run.1.log cytb.chiloglanis.run.2.log
logcombiner -trees -burnin 2000001 cytb.chiloglanis.run.1.trees cytb.chiloglanis.run.2.trees cytb.chiloglanis.run.1-2.trees
grep -c "tree STATE_" cytb.chiloglanis.run.1-2.trees
treeannotator -burninTrees 0 -heights ca cytb.chiloglanis.run.1-2.trees cytb.chiloglanis.beast.tre
```



```r
rm(list=ls())
today.dir <- glue('Results_{Sys.Date()}')
today.path <- here("temp",today.dir)
if(!dir.exists(today.path)) {dir.create(today.path,recursive=TRUE)}

# load data
raxml.tr.path <- here(today.path,"cytb.chiloglanis.rooted.nwk")
beast.tr.path <- here(today.path,"cytb.chiloglanis.rooted.tre")
cytb.chiloglanis.df <- read_csv(here("assets/cytb.chiloglanis.csv"),show_col_types=FALSE)
cytb.chiloglanis.fa.haps <- ape::read.FASTA(here(today.path,"cytb.chiloglanis.haps.fasta"))
cytb.chiloglanis.beast.tr <- ape::read.nexus(beast.tr.path)
cytb.chiloglanis.raxml.tr <- ape::read.tree(raxml.tr.path)

## run GMYC
# run gmyc simple
ape::is.binary(cytb.chiloglanis.beast.tr)
set.seed(42)
gmyc.res <- splits::gmyc(cytb.chiloglanis.beast.tr,method="single",interval=c(0,5),quiet=FALSE)
summary(gmyc.res)
# make df
gmyc.df <- gmyc_tbl(gmyc.res)
print(gmyc.df)

# run bgmyc
set.seed(42)
bgmyc.res.single <- bGMYC::bgmyc.singlephy(cytb.chiloglanis.beast.tr,mcmc=11000,burnin=1000,thinning=100,t1=2,t2=length(cytb.chiloglanis.beast.tr$tip.label),start=c(1,0.5,50))
# make df
bgmyc.df <- delimtools::bgmyc_tbl(bgmyc.res.single,ppcutoff=0.05)
print(bgmyc.df)

# run locmin
mat <- ape::dist.dna(cytb.chiloglanis.fa.haps,model="raw",pairwise.deletion=TRUE)
lmin <- spider::localMinima(as.matrix(mat))
locmin.df <- delimtools::locmin_tbl(mat,threshold=lmin$localMinima[1])
print(locmin.df)

#run asap
asap.df <- delimtools::asap(infile=here(today.path,"cytb.chiloglanis.haps.fasta"),model=3,outfolder=today.path)
print(asap.df)
# change groupe_1 to Partition_1

# run MPTP
# minbr
#minbrlen <- format(min(cytb.chiloglanis.raxml.tr.root$edge.length),scientific=FALSE)
#delimtools::minbr(tree=raxml.tr.path, file=here("assets/cytb.chiloglanis.fasta"))
mptp.df <- delimtools::mptp(infile=here(today.path,"cytb.chiloglanis.rooted.nwk"),outfolder=today.path,method="single")
print(mptp.df)
# fix locations of output files

# join
all.delims.df <- delim_join(list(gmyc.df,bgmyc.df,locmin.df,asap.df,mptp.df))

# autoplot
#cytb.chiloglanis.beast.tr.plot <- treeio::read.beast("assets/cytb.chiloglanis.rooted.tre")
cytb.chiloglanis.beast.tr.plot <- treeio::read.beast(beast.tr.path)


ftab <- cytb.chiloglanis.df |> 
    filter(gbAccession %in% pull(all.delims.df,labels)) |> 
    mutate(labs=glue("{gbAccession} | {scientificName}")) |> 
    dplyr::select(gbAccession,labs)

cols2 <- delim_brewer(all.delims.df,"Set1",9)


source(here("../delimtools/R/delim_autoplot2.R"))
library("ggtree")
cytb.chiloglanis.beast.tr.plot@data
p <- delim_autoplot2(delim=all.delims.df,tr=cytb.chiloglanis.beast.tr.plot,tbl_labs=ftab,col_vec=cols2,hexpand=1,widths=c(1,0.5)) 
ggsave(here(today.path,"delim.pdf"),plot=p,height=297,width=210,units="mm")

```