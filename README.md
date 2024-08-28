# Using delimtools to delimit Neotropical cichlids of the genus Geophagus

Please install required software following instructions in [install.md](assets/install.md).

For code to acquire the Geophagus dataset please follow instructions in [acquire-sequence-data.md](assets/acquire-sequence-data.md).


### R code to run Geophagus delimitation analysis using delimtools

```r
##################
#### packages ####
##################

#renv::install(here::here("../delimtools"))
#renv::install("legalLab/delimtools")
library("here")
library("glue")
library("tidyverse")
library("ape")
library("spider")
library("splits")
library("bGMYC")
library("treeio")
library("ggtree")
library("randomcoloR")
library("delimtools")


##################
###### setup #####
##################

# set temporary working directory in 'temp/'
rm(list=ls())
today.dir <- glue('Results_{Sys.Date()}')
today.path <- here("temp",today.dir)
if(!dir.exists(today.path)) {dir.create(today.path,recursive=TRUE)}

# fun to report n delims
report_delim <- function(df) {
    rep <- df |> pull(2) |> unique() |> length()
    df.name <- deparse(substitute(df))
    writeLines(glue::glue("The '{df.name}' delimitation table contains a total of {rep} unique species."))
}


##################
### load data ####
##################

coi.geophagus.haps.raxml.tr <- ape::read.tree(here("assets/coi.geophagus.haps.raxml.nwk"))
coi.geophagus.haps.beast.tr <- treeio::read.beast(here("assets/coi.geophagus.haps.beast.tre"))
coi.geophagus.haps.df <- readr::read_csv(here("assets/coi.geophagus.haps.csv"),show_col_types=FALSE)
coi.geophagus.haps.fa <- ape::read.FASTA(here("assets/coi.geophagus.haps.fasta"))


##################
#### run gmyc ####
##################

ape::is.binary(treeio::as.phylo(coi.geophagus.haps.beast.tr))
set.seed(42)
gmyc.res <- splits::gmyc(treeio::as.phylo(coi.geophagus.haps.beast.tr),method="single",interval=c(0,5),quiet=FALSE)
summary(gmyc.res)
# make df
gmyc.df <- gmyc_tbl(gmyc.res)
#gmyc.df |> print(n=Inf)
gmyc.df |> report_delim()


##################
### run bgmyc ###
##################

set.seed(42)
bgmyc.res.single <- bGMYC::bgmyc.singlephy(treeio::as.phylo(coi.geophagus.haps.beast.tr),mcmc=11000,burnin=1000,thinning=100,t1=2,t2=length(treeio::as.phylo(coi.geophagus.haps.beast.tr)$tip.label),start=c(1,0.5,50))
# make df
bgmyc.df <- delimtools::bgmyc_tbl(bgmyc.res.single,ppcutoff=0.05)
#bgmyc.df |> print(n=Inf)
bgmyc.df |> report_delim()


##################
### run locmin ###
##################

mat <- ape::dist.dna(coi.geophagus.haps.fa,model="raw",pairwise.deletion=TRUE)
lmin <- spider::localMinima(as.matrix(mat))
plot(lmin)
locmin.df <- delimtools::locmin_tbl(mat,threshold=lmin$localMinima[1])
#locmin.df |> print(n=Inf)
locmin.df |> report_delim()


##################
##### run 2% #####
##################

locmin.df.pc <- delimtools::locmin_tbl(mat,threshold=0.02) |> dplyr::rename(percent=locmin)
#locmin.df.pc |> print(n=Inf)
locmin.df.pc |> report_delim()


##################
#### run asap ####
##################

asap.df <- delimtools::asap(infile=here("assets/coi.geophagus.haps.fasta"),model=3,outfolder=today.path)
#asap.df |> print(n=Inf)
asap.df |> report_delim()


##################
#### run mptp ####
##################

#minbrlen <- format(min(coi.geophagus.raxml.tr.root$edge.length),scientific=FALSE)
#delimtools::minbr(tree=raxml.tr.path, file=here("assets/coi.geophagus.fasta"))
mptp.df <- delimtools::mptp(infile=here("assets/coi.geophagus.haps.raxml.nwk"),outfolder=today.path,method="single")
#mptp.df |> print(n=Inf)
mptp.df |> report_delim()


##################
### join delims ##
##################

all.delims.df <- delimtools::delim_join(list(gmyc.df,bgmyc.df,locmin.df,locmin.df.pc,asap.df,mptp.df))
#all.delims.df |> print(n=Inf)


##################
# subsample data #
##################

# clean and subsample
set.seed(42)
coi.geophagus.haps.df.sub <- coi.geophagus.haps.df |> 
    mutate(scientificName=str_replace_all(scientificName,"_AMX-2021","")) |> 
    mutate(scientificName=str_replace_all(scientificName,"_"," ")) |> 
    slice_sample(n=8,by=scientificName)

# subsample the delims
all.delims.df.sub <- all.delims.df |> filter(labels %in% pull(coi.geophagus.haps.df.sub,gbAccession))

# subample tips
#source(here("../delimtools/R/delim_autoplot2.R"))
coi.geophagus.haps.beast.tr.sub <- coi.geophagus.haps.beast.tr |> tidytree::keep.tip(pull(coi.geophagus.haps.df.sub,gbAccession))


##################
### plot delims ##
##################

# make tip label table
ftab <- coi.geophagus.haps.df.sub |> 
    mutate(labs=glue("{gbAccession} | {scientificName}")) |> 
    dplyr::select(gbAccession,labs)

# get number spp
#cols2 <- delim_brewer(all.delims.df,"Set1",n=9,seed=7)
n.spp <- all.delims.df.sub |> 
    pivot_longer(cols=!labels,names_to="method",values_to="spp") |> 
    distinct(spp) |> 
    pull(spp) |> 
    length()

# randomise colours
set.seed(42)
cols2 <- randomcoloR::distinctColorPalette(k=n.spp)

# plot and save
p <- delim_autoplot(delim=all.delims.df.sub,tr=coi.geophagus.haps.beast.tr.sub,tbl_labs=ftab,col_vec=cols2,hexpand=0.4,widths=c(0.5,0.2),n_match=3)
ggsave(here(today.path,"geophagus-delimitation.pdf"),plot=p,height=400,width=300,units="mm")
```
