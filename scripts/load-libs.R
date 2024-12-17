#!/usr/bin/env Rscript

# quiet messages
invisible(suppressMessages(suppressWarnings({
    library("here")
    library("glue")
    library("cli")
    library("knitr")
    library("tidyverse")
    library("ape")
    library("spider")
    library("splits")
    library("bGMYC")
    library("treeio")
    library("ggtree")
    library("randomcoloR")
    library("haplotypes")
    library("delimtools")
})))

pkk <- sessionInfo()
print(pkk)

# fun to report n delims
#report_delim <- function(df) {
#    rep <- df |> pull(2) |> unique() |> length()
#    df.name <- deparse(substitute(df))
#    writeLines(glue::glue("The '{df.name}' delimitation table contains a total of {rep} unique species."))
#}

 # report
writeLines("\n")
cli::cli_alert_success("{length(pkk$otherPkgs)} R packages loaded.")
writeLines("\n")
