#!/usr/bin/env Rscript

suppressMessages( library(tidyverse) )
library( optparse )
library( naivestates )

## Identify directory of the script
wd <- commandArgs( trailingOnly=FALSE ) %>%
    keep( ~grepl("--file=", .x) ) %>%
    str_replace( "--file=", "" ) %>% dirname()
cat( "Running the script from", wd, "\n" )

## Parse command-line arugments
option_list <- list(
    make_option(c("-i", "--in"), type="character", help="Input file"),
    make_option(c("-o", "--out"), type="character", default="/data",
                help="Output directory"),
    make_option(c("-m", "--markers"), type="character", default="auto",
                help="Markers to model"),
    make_option(c("-p", "--plots"), type="character", default="off",
                help="Generate plots showing the fit"),
    make_option("--mct", type="character", default="",
                help="Marker -> cell type map in .csv format"),
    make_option("--id", type="character", default="CellID",
                help="Column containing cell IDs"),
    make_option("--log", type="character", default="auto",
                help="Whether to apply a log transform <yes|no|auto>"),
    make_option("--sfx", type="character", default="",
                help="Common suffix on marker columns (e.g., _cellMask)"),
    make_option("--umap", action="store_true", default=FALSE,
                help="Generate UMAP plots")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Argument verification
if( !("in" %in% names(opt)) )
    stop( "Please provide an input file name with -i" )
if( !(opt$log %in% c("yes","no","auto")) )
    stop( "--log must be one of <yes|no|auto>" )
if( !(opt$plots %in% c("off", "pdf", "png")) )
    stop( "--plots must be one of <off|pdf|png>" )

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

## Identify markers in the matrix
mrkv <- findMarkers(setdiff(colnames(X), opt$id), opt$markers,
                    opt$sfx, TRUE, TRUE, TRUE)

## Handle log transformation of the data
if( opt$log == "yes" ||
    (opt$log == "auto" && max(X[mrkv], na.rm=TRUE) > 1000) )
{
    cat( "Applying a log10 transform\n" )
    X <- X %>% mutate_at( unname(mrkv), ~log10(.x+1) )
}

## Fit Gaussian mixture models
GMM <- GMMfit(X, opt$id, !!!mrkv)
fnMdl <- file.path( opt$out, str_c(sn, "-models.csv") )
cat( "Saving models to", fnMdl, "\n" )
GMMmodels(GMM) %>% write_csv( fnMdl )

## Reshape the matrix back to cells-by-marker format
Y <- GMMreshape(GMM)

cat( "------\n" )

## Find the default cell type map
if( opt$mct != "" ) {

    ## Load marker -> cell type associations
    cat( "Loading cell type map from", opt$mct, "\n" )
    mct <- read_csv( opt$mct, col_types=cols() ) %>%
        distinct() %>% filter(Marker %in% colnames(Y))

    if( nrow(mct) == 0 ) {
        warning( "No usable marker -> cell type mappings detected" )
        Y <- findDominant(Y, opt$id)
    } else {
        cat( "Using the following marker -> cell type map:\n" )
        walk2( mct$Marker, mct$State, ~cat(.x, "->", .y, "\n") )
        Y <- callStates(Y, opt$id, mct)
    }
} else {
    cat( "No marker -> cell type mapping provided\n" )
    Y <- findDominant(Y, opt$id)
}

cat( "------\n" )

## Identify the output location(s)
fnOut <- file.path( opt$out, str_c(sn, "-states.csv") )
cat( "Saving probabilities and calls to", fnOut, "\n")
Y %>% write_csv( fnOut )

## Generates plots as necessary
if( opt$plots != "off" )
{
    ## Create a separate directory for plots
    dirPlot <- file.path( opt$out, "plots", sn )
    dir.create(dirPlot, recursive=TRUE, showWarnings=FALSE)

    ## Fit overview
    fn <- file.path( file.path(opt$out, "plots"), str_c(sn, "-allfits.", opt$plots) )
    ggf <- plotFitOverview(GMM)
    suppressMessages(ggsave( fn, ggf, width=12, height=8 ))
    
    ## Compute a UMAP projection
    if( opt$umap ) {
        cat( "Computing a UMAP projection...\n" )
        U <- umap( Y, c(opt$id, "State", "Dominant") )
    
        ## Generate and write a summary plot
        gg <- plotSummary( U )
        fn <- file.path( file.path(opt$out, "plots"), str_c(sn, "-summary.", opt$plots) )
        suppressMessages(ggsave( fn, gg, width=9, height=7 ))
        cat( "Plotted summary to", fn, "\n" )

        ## Generate and write faceted probabilities plot
        gg <- plotProbs( U, c(opt$id, "State", "Dominant") )
        fn <- file.path( file.path(opt$out, "plots"), str_c(sn, "-probs.", opt$plots) )
        suppressMessages(ggsave( fn, gg, width=9, height=7 ))
        cat( "Plotted probabilities to", fn, "\n" )
    }

    ## Generate and write out plots for individual marker fits
    for( i in names(mrkv) )
    {
        gg <- plotMarker(GMM, i)
        fn <- file.path( dirPlot, str_c(i,".",opt$plots) )
        suppressMessages(ggsave( fn, gg ))
        cat( "Wrote", fn, "\n" )
    }
}
