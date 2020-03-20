#!/usr/bin/env Rscript

suppressMessages( library(tidyverse) )
library( optparse )
library( naivestates )

## Parse command-line arugments
option_list <- list(
    make_option(c("-i", "--in"), type="character", help="Input file"),
    make_option(c("-o", "--out"), type="character", default=".",
                help="Output directory"),
    make_option(c("-m", "--markers"), type="character", default="DNA0",
                help="Markers to model"),
    make_option(c("-p", "--plots"), action="store_true", default=FALSE,
                help="Generate plots showing the fit"),
    make_option(c("--id"), type="character", default="CellId",
                help="Column containing cell IDs")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Argument verification
if( !("in" %in% names(opt)) )
    stop( "Please provide an input file name with -i" )

## Identify the sample name
sn <- basename( opt[["in"]] ) %>% str_split( "\\." ) %>%
    pluck( 1, 1 )
cat( "Inferred sample name:", sn, "\n" )

## Read the data matrix
X <- read_csv( opt[["in"]], col_types=cols() )
cat( "Read", nrow(X), "entries\n" )

## Identify markers in the matrix
mrk <- str_split( opt$markers, "," )[[1]] %>% set_names()
mrki <- map( mrk, grep, colnames(X) )

## Verify marker uniqueness
iwalk( mrki, ~if(length(.x) > 1)
                  stop("Marker ", .y, " maps to multiple columns") )
iwalk( mrki, ~if(length(.x) == 0)
                  stop("Marker ", .y, " is not found in the data") )

## Fit Gaussian mixture models
GMM <- GMMfit( X, opt$id, !!!mrki, baseline=0.01 )

## Identify the output location(s)
fnOut <- file.path( opt$out, str_c(sn, "_ep.csv") )
cat( "Saving expression probabilities to", fnOut, "\n")
GMMreshape(GMM) %>% write_csv( fnOut )

## Generates plots as necessary
if( opt$plots )
{
    ## Create a separate directory for plots
    dirPlot <- file.path( opt$out, "plots", sn )
    dir.create(dirPlot, recursive=TRUE, showWarnings=FALSE)

    ## Generate and write out individual plots
    for( i in names(mrki) )
    {
        gg <- plotFit(GMM, i)
        fn <- file.path( dirPlot, str_c(i,".pdf") )
        ggsave( fn, gg )
        cat( "Wrote", fn, "\n" )
    }
}
