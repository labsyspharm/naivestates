#!/usr/bin/env Rscript

suppressMessages( library(tidyverse) )
library( optparse )
library( naivestates )

## Parse command-line arugments
option_list <- list(
    make_option(c("-i", "--in"), type="character", help="Input file"),
    make_option(c("-o", "--out"), type="character", default="/data",
                help="Output directory"),
    make_option(c("-m", "--markers"), type="character", default="auto",
                help="Markers to model"),
    make_option(c("-p", "--plots"), action="store_true", default=FALSE,
                help="Generate plots showing the fit"),
    make_option("--mct", type="character", default="/app/typemap.csv",
                help="Marker -> cell type map in .csv format"),
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
if( file.exists(opt$markers) ) {
    mrk <- scan(opt$markers, what=character(), quiet=TRUE)
} else if( opt$markers == "auto" ) {
    mrk <- autoMarkers(setdiff(colnames(X), opt$id))
} else {
    mrk <- str_split( opt$markers, "," )[[1]]
}

## Identify markers in the matrix
mrkv <- findMarkers( colnames(X), mrk, TRUE, TRUE )
cat( "Found markers:", str_flatten(names(mrkv), ", "), "\n" )

## Handle log transformation of the data
if( opt$log == "yes" ||
    (opt$log == "auto" && max(X[mrkv]) > 1000) )
{
    cat( "Applying a log10 transform\n" )
    X <- X %>% mutate_at( mrkv, ~log10(.x+1) )
}

## Fit Gaussian mixture models
GMM <- GMMfit(X, opt$id, !!!mrkv)
Y <- GMMreshape(GMM)

cat( "------\n" )

## Load marker -> cell type associations
tm <- read_csv( opt$mct, col_types=cols() ) %>% deframe()
mct <- findMarkers( colnames(Y), names(tm) )
mct <- set_names( tm[names(mct)], mct )

if( length(mct) == 0 ) {
    warning( "No usable marker -> cell type mappings detected" )
    Y <- callStates(Y, opt$id)
} else {
    cat( "Using the following marker -> cell type map:\n" )
    iwalk( mct, ~cat( .y, "->", .x, "\n" ) )
    Y <- callStates(Y, opt$id, mct)
}

cat( "------\n" )

## Identify the output location(s)
fnOut <- file.path( opt$out, str_c(sn, "-states.csv") )
cat( "Saving results to", fnOut, "\n")
Y %>% write_csv( fnOut )

## Generates plots as necessary
if( opt$plots )
{
    ## Create a separate directory for plots
    dirPlot <- file.path( opt$out, "plots", sn )
    dir.create(dirPlot, recursive=TRUE, showWarnings=FALSE)

    ## Compute a UMAP projection
    U <- umap( Y, c(opt$id, "State", "Anchor") )
    
    ## Generate and write a summary plot
    gg <- plotSummary( U )
    fn <- file.path( file.path(opt$out, "plots"), str_c(sn, "-summary.pdf") )
    suppressMessages(ggsave( fn, gg, width=9, height=7 ))
    cat( "Wrote summary to", fn, "\n" )

    ## Generate and write out plots for individual marker fits
    for( i in names(mrkv) )
    {
        gg <- plotFit(GMM, i)
        fn <- file.path( dirPlot, str_c(i,".pdf") )
        suppressMessages(ggsave( fn, gg ))
        cat( "Wrote", fn, "\n" )
    }
}
