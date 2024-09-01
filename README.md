# Using delimtools to delimit Neotropical cichlids of the genus _Geophagus_

Please install required software following instructions in [install.md](scripts/install.md).

For code to acquire the _Geophagus_ dataset please follow instructions in [acquire-sequence-data.md](scripts/acquire-sequence-data.md).


### R code to run _Geophagus_ delimitation analysis using delimtools

```r
##################
#### packages ####
##################

#renv::install(here::here("../delimtools"))
#renv::install("legalLab/delimtools")
# load R packages
source(here::here("scripts/load-libs.R"))


##################
###### setup #####
##################

# set temporary working directory in 'temp/'
today.dir <- glue::glue('Results_{Sys.Date()}')
today.path <- here::here("temp",today.dir)
if(!dir.exists(today.path)) {dir.create(today.path,recursive=TRUE)}



##################
### load data ####
##################

coi.geophagus.haps.raxml.tr <- ape::read.tree(here::here("assets/coi.geophagus.haps.raxml.nwk"))
coi.geophagus.haps.beast.tr <- treeio::read.beast(here::here("assets/coi.geophagus.haps.beast.tre"))
coi.geophagus.haps.df <- readr::read_csv(here::here("assets/coi.geophagus.haps.csv"),show_col_types=FALSE)
coi.geophagus.haps.fa <- ape::read.FASTA(here::here("assets/coi.geophagus.haps.fasta"))


##################
#### run gmyc ####
##################

# check it tree is binary (should be TRUE)
ape::is.binary(treeio::as.phylo(coi.geophagus.haps.beast.tr))
set.seed(42)
gmyc.res <- splits::gmyc(treeio::as.phylo(coi.geophagus.haps.beast.tr),method="single",interval=c(0,5),quiet=FALSE)
summary(gmyc.res)
# make df
gmyc.df <- delimtools::gmyc_tbl(gmyc.res)
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
plot(lmin); abline(v=lmin$localMinima[1],col="red")
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

source(here("../delimtools/R/asap.R"))
asap.df <- asap(infile=here::here("assets/coi.geophagus.haps.fasta"),exe=here::here("software/ASAP/bin/asap"),model=3)
#asap.df <- delimtools::asap(infile=here::here("assets/coi.geophagus.haps.fasta"),exe=here::here("software/ASAP/bin/asap"),model=3)
#asap.df |> print(n=Inf)
asap.df |> report_delim()


##################
### run mptp s ###
##################

#minbrlen <- format(min(coi.geophagus.raxml.tr.root$edge.length),scientific=FALSE)
#delimtools::minbr(tree=raxml.tr.path, file=here("assets/coi.geophagus.fasta"))
#source(here("../delimtools/R/mptp.R"))
#mptp.s.df <- mptp(infile=here("assets/coi.geophagus.haps.raxml.nwk"),exe=here::here("software/mptp/bin/mptp"),method="single")
mptp.s.df <- delimtools::mptp(infile=here("assets/coi.geophagus.haps.raxml.nwk"),exe=here::here("software/mptp/bin/mptp"),method="single")
#mptp.df |> print(n=Inf)
mptp.s.df |> report_delim()


##################
### run mptp m ###
##################

#minbrlen <- format(min(coi.geophagus.raxml.tr.root$edge.length),scientific=FALSE)
#delimtools::minbr(tree=raxml.tr.path, file=here("assets/coi.geophagus.fasta"))
mptp.m.df <- delimtools::mptp(infile=here("assets/coi.geophagus.haps.raxml.nwk"),exe=here::here("software/mptp/bin/mptp"),outfolder=today.path,method="multi")
#mptp.df |> print(n=Inf)
mptp.m.df |> report_delim()


##################
### join delims ##
##################

all.delims.df <- delimtools::delim_join(list(gmyc.df,bgmyc.df,locmin.df,locmin.df.pc,asap.df,mptp.s.df,mptp.m.df))
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
tip.tab <- coi.geophagus.haps.df.sub |> 
    dplyr::mutate(labs=glue::glue("{gbAccession} | {scientificName}")) |> 
    dplyr::select(gbAccession,labs)

# get number spp
#cols2 <- delim_brewer(all.delims.df,"Set1",n=9,seed=7)
n.spp <- all.delims.df.sub |> 
    tidyr::pivot_longer(cols=!labels,names_to="method",values_to="spp") |> 
    dplyr::distinct(spp) |> 
    dplyr::pull(spp) |> 
    length()

# randomise colours
set.seed(42)
cols2 <- randomcoloR::distinctColorPalette(k=n.spp)

# plot and save
p <- delimtools::delim_autoplot(delim=all.delims.df.sub,tr=coi.geophagus.haps.beast.tr.sub,tbl_labs=ftab,col_vec=cols2,hexpand=0.4,widths=c(0.5,0.2),n_match=3,delim_order=c("asap","locmin","percent","gmyc","bgmyc","ptp","mptp"))
ggplot2::ggsave(here::here(today.path,"geophagus-delimitation.pdf"),plot=p,height=400,width=300,units="mm")
```
