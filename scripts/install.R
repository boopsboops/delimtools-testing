#!/usr/bin/env Rscript

# install renv
install.packages("renv")
library("renv")

# install delimtools bioconductor dependencies
renv::install(c("bioc::BiocVersion","bioc::ggtree","bioc::ggtreeExtra","bioc::treeio"))

# install local packages
renv::install("here")
library("here")
renv::install(here("../splits_1.0-19.tar.gz"))
renv::install(here("../bGMYC_1.0.2.tar.gz"))

# install delimtools
# need to check github PAT is available
# see https://usethis.r-lib.org/articles/git-credentials.html
gitcreds::gitcreds_set()
usethis::gh_token_help()
renv::install("legalLab/delimtools")
# can specify a specific commit 
renv::install("legalLab/delimtools@c94e41a")

# record packages
renv::status()
renv::snapshot()
