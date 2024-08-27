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
scripts/download-sequences.R -c Geophagus -n 500 -x 2500 -b 1 -a false -d true
# real
scripts/download-sequences.R -c Geophagus -n 500 -x 2500 -b 500 -a false -d false
# cluster
scripts/clean-and-cluster.R -n 10 -c 0.6 -m 2
# pick
scripts/pick-clusters.R -c 10 -g coi
# annotate
scripts/annotate-ncbi.R -t 3 -c fishbase
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

# rm(list=ls())
today.dir <- glue('Results_{Sys.Date()}')
today.path <- here("temp",today.dir)
if(!dir.exists(today.path)) {dir.create(today.path,recursive=TRUE)}

# # copy
file.copy(here("ncbi-supermatrix/temp/Results_2024-08-27/alignments/coi.aligned.trimmed.fasta"), here(today.path,"coi.geophagus.fasta"))
file.copy(here("ncbi-supermatrix/temp/Results_2024-08-27/trees/coi.aligned.trimmed.fasta.raxml.bestTree"), here(today.path,"coi.geophagus.nwk"))
file.copy(here("ncbi-supermatrix/temp/Results_2024-08-27/metadata/ncbi-clean.csv"), here(today.path,"coi.geophagus.csv"))

# read fasta
coi.geophagus.fa <- ape::read.FASTA(here(today.path,"coi.geophagus.fasta"))

# read tree
coi.geophagus.raxml.tr <- ape::read.tree(here(today.path,"coi.geophagus.nwk"))

# read data 
coi.geophagus.df <- read_csv(here(today.path,"coi.geophagus.csv"),show_col_types=FALSE)

# get ingroup and outgroups
ingroup.tips <- coi.geophagus.df |> 
    filter(grepl("harreri|winemilleri",scientificName)) |> 
    slice_head(n=1,by=scientificName) |> 
    pull(gbAccession)

# get one outgroup
outgroup <- coi.geophagus.df |> 
    filter(grepl("brasiliensis",scientificName)) |> 
    slice_head(n=1) |> 
    pull(gbAccession)

# get all ingroup tips
coi.geophagus.raxml.tr.rooted <- coi.geophagus.raxml.tr |> ape::root(outgroup=outgroup,resolve.root=TRUE)
ingroup.node <- coi.geophagus.raxml.tr.rooted |> ape::getMRCA(ingroup.tips)
ingroup <- extract.clade(coi.geophagus.raxml.tr.rooted,node=ingroup.node)$tip.label

# join ingroup and outgroup
geo.gb <- c(ingroup,outgroup)

#coi.geophagus.beast.tr <- ape::read.nexus(here("temp/beast/coi.geophagus.beast.tre"))
#coi.geophagus.df |> glimpse()

# subset ingroup from fasta
coi.geophagus.ingroup.fa <- coi.geophagus.fa[geo.gb]

# collapse haplotypes
# hap table
delimtools::haplotype_tbl(coi.geophagus.ingroup.fa,verbose=FALSE) |> print(n=Inf)
haps.collapsed <- delimtools::hap_collapse(coi.geophagus.ingroup.fa,collapseSubstrings=TRUE,clean=TRUE) |> names()

# subset and write out
coi.geophagus.haps.fa <- coi.geophagus.ingroup.fa[haps.collapsed]
coi.geophagus.haps.fa |> ape::write.FASTA(here(today.path,"coi.geophagus.haps.fasta"))

# write out df
coi.geophagus.df |> 
    filter(gbAccession %in% haps.collapsed) |> 
    write_csv(here(today.path,"coi.geophagus.haps.csv"))
```


```bash
# get the current temp dir
tmpdir=$(ls -d temp/*/ | sort -r | head -n 1)
echo $tmpdir
cd $tmpdir
beauti coi.geophagus.haps.fasta

# run beast
beast -beagle_auto -overwrite -seed 42 coi.geophagus.haps.run1.xml 
beast -beagle_auto -overwrite -seed 24 coi.geophagus.haps.run2.xml 
tracer coi.geophagus.haps.run1.log coi.geophagus.haps.run2.log
logcombiner -trees -burnin 2016000 coi.geophagus.haps.run1.trees coi.geophagus.haps.run2.trees coi.geophagus.haps.run1+2.trees
grep -c "tree STATE_" coi.geophagus.haps.run1+2.trees
treeannotator -burninTrees 0 -heights ca coi.geophagus.haps.run1+2.trees coi.geophagus.haps.beast.tre

# run raxml on beast tree
# run in R to convert tree
#ape::read.nexus(here(today.path,"coi.geophagus.haps.beast.tre")) |> write.tree(here(today.path,"coi.geophagus.haps.beast.tre.nwk"))
raxml-ng --evaluate --threads auto{} --tree coi.geophagus.haps.beast.tre.nwk --lh-epsilon 0.1 --redo --seed 42 --outgroup MH538063.1 --model TN93+G --msa coi.geophagus.haps.fasta
```



```r
# get temp dir
today.dir <- glue('Results_{Sys.Date()}')
today.path <- here("temp",today.dir)
if(!dir.exists(today.path)) {dir.create(today.path,recursive=TRUE)}

# read tree and drop outgroup
coi.geophagus.haps.raxml.tr <- ape::read.tree(here(today.path,"coi.geophagus.haps.fasta.raxml.bestTree")) |> ape::drop.tip("MH538063.1")

# read beast tree and drop outgroup
coi.geophagus.haps.beast.tr <- treeio::read.beast(here(today.path,"coi.geophagus.haps.beast.tre")) |> treeio::drop.tip("MH538063.1")
#treeio::as.phylo(coi.geophagus.haps.beast.tr)

# read data 
coi.geophagus.haps.df <- read_csv(here(today.path,"coi.geophagus.haps.csv"),show_col_types=FALSE) |> filter(gbAccession!="MH538063.1")

# read dna 
coi.geophagus.haps.fa <- ape::read.FASTA(here(today.path,"coi.geophagus.haps.fasta"))
coi.geophagus.haps.fa <- coi.geophagus.haps.fa[which(names(coi.geophagus.haps.fa)!="MH538063.1")]

# write out
coi.geophagus.haps.raxml.tr |> ape::write.tree(here("assets/coi.geophagus.haps.raxml.nwk"))
coi.geophagus.haps.beast.tr |> treeio::write.beast(here("assets/coi.geophagus.haps.beast.tre"))
coi.geophagus.haps.df |> readr::write_csv(here("assets/coi.geophagus.haps.csv"))
coi.geophagus.haps.fa |> ape::write.FASTA(here("assets/coi.geophagus.haps.fasta"))

######################
######################
```


```r
# set dirs
rm(list=ls())
today.dir <- glue('Results_{Sys.Date()}')
today.path <- here("temp",today.dir)
if(!dir.exists(today.path)) {dir.create(today.path,recursive=TRUE)}

# fun to report
report_delim <- function(df) {
    rep <- df |> pull(2) |> unique() |> length()
    df.name <- deparse(substitute(df))
    writeLines(glue::glue("The '{df.name}' delimitation table contains a total of {rep} unique species."))
}


# read in
coi.geophagus.haps.raxml.tr <- ape::read.tree(here("assets/coi.geophagus.haps.raxml.nwk"))
coi.geophagus.haps.beast.tr <- treeio::read.beast(here("assets/coi.geophagus.haps.beast.tre"))
coi.geophagus.haps.df <- readr::read_csv(here("assets/coi.geophagus.haps.csv"),show_col_types=FALSE)
coi.geophagus.haps.fa <- ape::read.FASTA(here("assets/coi.geophagus.haps.fasta"))

## run GMYC
# run gmyc simple
ape::is.binary(treeio::as.phylo(coi.geophagus.haps.beast.tr))
set.seed(42)
gmyc.res <- splits::gmyc(treeio::as.phylo(coi.geophagus.haps.beast.tr),method="single",interval=c(0,5),quiet=FALSE)
summary(gmyc.res)
# make df
gmyc.df <- gmyc_tbl(gmyc.res)
print(n=Inf)
gmyc.df |> report_delim()


# run bgmyc
set.seed(42)
bgmyc.res.single <- bGMYC::bgmyc.singlephy(treeio::as.phylo(coi.geophagus.haps.beast.tr),mcmc=11000,burnin=1000,thinning=100,t1=2,t2=length(treeio::as.phylo(coi.geophagus.haps.beast.tr)$tip.label),start=c(1,0.5,50))
# make df
bgmyc.df <- delimtools::bgmyc_tbl(bgmyc.res.single,ppcutoff=0.05)
bgmyc.df |> print(n=Inf)
bgmyc.df |> report_delim()

# run locmin
mat <- ape::dist.dna(coi.geophagus.haps.fa,model="raw",pairwise.deletion=TRUE)
lmin <- spider::localMinima(as.matrix(mat))
plot(lmin)
locmin.df <- delimtools::locmin_tbl(mat,threshold=lmin$localMinima[1])
locmin.df |> print(n=Inf)
locmin.df |> report_delim()

# run locmin 2%
locmin.df.pc <- delimtools::locmin_tbl(mat,threshold=0.02) |> dplyr::rename(percent=locmin)
locmin.df.pc |> print(n=Inf)
locmin.df.pc |> report_delim()

#run asap
asap.df <- delimtools::asap(infile=here("assets/coi.geophagus.haps.fasta"),model=3,outfolder=today.path)
asap.df |> print(n=Inf)
asap.df |> report_delim()


# run MPTP
# minbr
#minbrlen <- format(min(coi.geophagus.raxml.tr.root$edge.length),scientific=FALSE)
#delimtools::minbr(tree=raxml.tr.path, file=here("assets/coi.geophagus.fasta"))
mptp.df <- delimtools::mptp(infile=here("assets/coi.geophagus.haps.raxml.nwk"),outfolder=today.path,method="single")
print(mptp.df,n=Inf)
mptp.df |> report_delim()

# join
all.delims.df <- delimtools::delim_join(list(gmyc.df,bgmyc.df,locmin.df,locmin.df.pc,asap.df,mptp.df))
all.delims.df |> print(n=Inf)

# subset
set.seed(42)
coi.geophagus.haps.df.sub <- coi.geophagus.haps.df |> 
    mutate(scientificName=str_replace_all(scientificName,"_AMX-2021","")) |> 
    mutate(scientificName=str_replace_all(scientificName,"_"," ")) |> 
    slice_sample(n=8,by=scientificName)

# subsample the delims
all.delims.df.sub <- all.delims.df |> filter(labels %in% pull(coi.geophagus.haps.df.sub,gbAccession))

# 
ftab <- coi.geophagus.haps.df.sub |> 
    #filter(gbAccession %in% pull(all.delims.df,labels)) |> 
    mutate(labs=glue("{gbAccession} | {scientificName}")) |> 
    dplyr::select(gbAccession,labs)

# choose cols
#cols2 <- delim_brewer(all.delims.df,"Set1",n=9,seed=7)

n.spp <- all.delims.df.sub |> pivot_longer(cols=!labels,names_to="method",values_to="spp") |> distinct(spp) |> pull(spp) |> length()
set.seed(42)
cols2 <- randomcoloR::distinctColorPalette(k=n.spp)


#source(here("../delimtools/R/delim_autoplot2.R"))
#library("ggtree")
coi.geophagus.haps.beast.tr.sub <- coi.geophagus.haps.beast.tr |> tidytree::keep.tip(pull(coi.geophagus.haps.df.sub,gbAccession))

p <- delim_autoplot(delim=all.delims.df.sub,tr=coi.geophagus.haps.beast.tr.sub,tbl_labs=ftab,col_vec=cols2,hexpand=0.4,widths=c(0.5,0.2),n_match=3) #widths=c(1,0.3)
ggsave(here(today.path,"delim.pdf"),plot=p,height=400,width=300,units="mm")

```
