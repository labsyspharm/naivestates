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
    make_option("--id", type="character", default="CellID",
                help="Column containing cell IDs"),
    make_option("--log", type="character", default="auto",
                help="Whether to apply a log transform <yes|no|auto>")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Argument verification
if( !("in" %in% names(opt)) )
    stop( "Please provide an input file name with -i" )
if( !(opt$log %in% c("yes","no","auto")) )
    stop( "--log must be one of <yes|no|auto>" )

## Identify the sample name
sn <- basename( opt$`in` ) %>% str_split( "\\." ) %>%
    pluck( 1, 1 )
cat( "Inferred sample name:", sn, "\n" )

## Read the data matrix
X <- read_csv( opt$`in`, col_types=cols() )
cat( "Read", nrow(X), "entries\n" )

## Fix potential capitalization mismatch of --id
if( !(opt$id %in% colnames(X)) )
{
    ## Attempt to find a singular case-insensitive match
    i <- grep( tolower(opt$id), tolower(colnames(X)) )
    if( length(i) == 1 )
    {
        warning( "  No such column ", opt$id,
                "; using ", colnames(X)[i], " instead" )
        opt$id <- colnames(X)[i]
    }
    else stop( "No such column ", opt$id,
              "; use --id to specify which column contains cell IDs" )
}

## Determine if we're working with a file of markers or if
##   markers are specified as a comma,delimited,list
mrk <- `if`( file.exists(opt$markers),
            scan(opt$markers, what=character(), quiet=TRUE),
            str_split( opt$markers, "," )[[1]] ) %>%
    set_names() %>% map_chr( str_c, "$" )

## Identify markers in the matrix
mrki <- map( mrk, grep, colnames(X) )
mrkv <- unlist(mrki)

## Verify marker uniqueness
iwalk( mrki, ~if(length(.x) > 1)
                  stop("Marker ", .y, " maps to multiple columns") )
iwalk( mrki, ~if(length(.x) == 0)
                  stop("Marker ", .y, " is not found in the data") )
cat( "Found markers:", str_flatten(names(mrki), ", "), "\n" )

## Handle log transformation of the data
if( opt$log == "yes" ||
    (opt$log == "auto" && max(X[mrkv]) > 1000) )
{
    cat( "Applying a log10 transform\n" )
    X <- X %>% mutate_at( colnames(X)[mrkv], ~log10(.x+1) )
}

## Fit Gaussian mixture models
GMM <- GMMfit(X, opt$id, !!!mrki)

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
