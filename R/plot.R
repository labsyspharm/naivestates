## Plotting functions
##
## by Artem Sokolov

## Bold element_text of desired font size s
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
          legend.title = etxt(14), strip.text = etxt(12),
          strip.background = element_blank() )
}

## Generates a palette for a vector cell type labels
##
## Not exported
makePal <- function( v )
{
    ## Visually distinct 20 colors + black
    ## by Sasha Trubetskoy
    vd20 <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8",
              "#f58231", "#911eb4", "#42d4f4", "#f032e6",
              "#bfef45", "#fabed4", "#469990", "#dcbeff",
              "#9A6324", "#fffac8", "#800000", "#aaffc3",
              "#808000", "#ffd8b1", "#000075", "#a9a9a9",
              "#000000")
    
    v1 <- setdiff( v, "None" ) %>% sort
    rlang::set_names( vd20[1:length(v1)], v1 ) %>%
        c( "None" = "lightgray" )
}

#' Plots a GMM fit for a single marker
#'
#' @param FT - the data frame returned by GMMfit()
#' @param marker - name of the marker to examine
#' @param plotFit - whether to plot the fit
#' @param plotMap - whether to plot expression -> probability map
#' @return a grid object containing a two-panel plot
#' @export
plotMarker <- function( FT, marker, plotFit=TRUE, plotMap=TRUE )
{
    ## Identify the slice to plot
    X1 <- FT %>% dplyr::filter( Marker == marker ) %>% purrr::pluck("Values",1) %>%
        na.omit() %>% dplyr::filter( AdjVal >= 0, AdjVal <= 1 ) %>% dplyr::arrange(Value)

    ## Identify the matching GMM and the plotting range
    gmm <- FT %>% dplyr::filter( Marker == marker ) %>% purrr::pluck("GMM", 1)
    xrng <- seq(gmm$lo, gmm$hi, length.out=5)
    xtck <- round(xrng, 2)

    ## Density plot
    gg1 <- ggplot2::ggplot( X1, aes(x=AdjVal) ) + ggplot2::theme_bw() +
        bold_theme() +
        ggplot2::ylab("Density") +
        ggplot2::xlab("Log-transformed Expression") +
        ggplot2::geom_density(lwd=1.1) +
        ggplot2::facet_wrap( ~str_c(marker) ) +
        ggplot2::scale_x_continuous(breaks = seq(0,1,by=0.25),
                                    labels = xtck)
    if( plotFit )
        gg1 <- gg1 +
            ggplot2::geom_line( aes(y=CP), color="red", lwd=1.1 ) +
            ggplot2::geom_line( aes(y=CN), color="blue", lwd=1.1 )
    
    if( !plotMap ) return(gg1)

    ## Remove duplicated axis information
    ebl <- ggplot2::element_blank()
    gg1 <- gg1 +
        ggplot2::theme(axis.title.x = ebl,
                       axis.text.x  = ebl,
                       axis.ticks.x = ebl)

    ## Probability plot
    gg2 <- ggplot2::ggplot( X1, aes(x=Value, y=Prob) ) +
        ggplot2::geom_line(lwd=1.1) + ggplot2::theme_bw() +
        bold_theme() +
        ggplot2::ylab( "P( Expressed )" ) + ggplot2::xlab("Log-transformed Expression") +
        ggplot2::scale_y_continuous( breaks = seq(0,1,by=0.2), limits=c(0,1) ) +
        ggplot2::scale_x_continuous( breaks = xrng, labels = xtck )
    egg::ggarrange( gg1, gg2, nrow=2 )
}

#' Plots an overview of all fits
#'
#' @param FT - the data frame returned by GMMfit()
#' @param npt - downsample to this many points (default: 100,000)
#' @importFrom magrittr %>%
#' @return a faceted ggplot object showing a fit to every marker
#' @export
plotFitOverview <- function( FT, npt=100000 )
{
    ## Downsample only if the number of points exceeds the request
    vs <- 1:nrow(FT$Values[[1]])
    if( length(vs) > npt ) vs <- sample(vs, npt)

    X <- FT %>% dplyr::mutate_at( "Values", map, slice, vs ) %>%
        dplyr::select( Marker, Values ) %>% tidyr::unnest( Values ) %>%
        dplyr::filter( AdjVal >= 0, AdjVal <= 1 )
    ggplot2::ggplot( X, ggplot2::aes(x=AdjVal) ) +
        ggplot2::theme_bw() + bold_theme() +
            ggplot2::ylab( "Density" ) +
            ggplot2::xlab( "Normalized Expression" ) +
            ggplot2::geom_density(lwd=1.1) +
            ggplot2::facet_wrap( ~Marker, scales="free_y" ) + 
            ggplot2::geom_line( ggplot2::aes(y=CP), color="red", lwd=1.1 ) +
            ggplot2::geom_line( ggplot2::aes(y=CN), color="blue", lwd=1.1 )
}

#' Summary UMAP plot with points colored by cell type (dominant label)
#' @param U data frame produced by umap()
#' @return ggplot object encapsulating the plot
#' @export
plotSummary <- function( U )
{
    v <- c("State", "Anchor") %in% colnames(U)
    if( "State" %in% colnames(U) ) {
        ggplot2::ggplot( U, ggplot2::aes(UMAP1, UMAP2, color=State) ) +
            ggplot2::geom_point() + ggplot2::theme_bw() +
                ggplot2::scale_color_manual(values=makePal(U$State),
                                            name="Cell Type/State")
    } else if( ("Dominant" %in% colnames(U)) &&
               length(unique(U$Dominant)) <= 21 ) {
        ## Plot UMAP projection, coloring by cell state calls
        ggplot2::ggplot( U, ggplot2::aes(UMAP1, UMAP2, color=Dominant) ) +
            ggplot2::geom_point() + ggplot2::theme_bw() +
                ggplot2::scale_color_manual(values=makePal(U$Dominant),
                                            name="Dominant Marker" )
    } else {
        ggplot2::ggplot( U, aes(UMAP1, UMAP2) ) +
            ggplot2::theme_bw() + ggplot2::geom_point(color="gray")
    }
}

#' Faceted UMAP plot showing probabilities of individual markers
#' @param U data frame produced by umap()
#' @param excl character vector of column names to exclude from plotting
#' @return ggplot object of the faceted plot
#' @importFrom magrittr %>%
#' @export
plotProbs <- function( U, excl )
{
    trm <- intersect( excl, colnames(U) )
    Z <- U %>% dplyr::select( -dplyr::one_of(trm) ) %>%
        tidyr::gather( Marker, Probability, -UMAP1, -UMAP2 )

    nmrk <- length( unique(Z$Marker) )
    
    ggplot2::ggplot( Z, ggplot2::aes(UMAP1, UMAP2, color=Probability) ) +
        ggplot2::geom_point(size=0.5) + ggplot2::theme_void() +
        ggplot2::facet_wrap( ~Marker, ncol=round(sqrt(nmrk)) ) +
        ggplot2::scale_color_gradient2(midpoint=0.5,
                                       limits=c(0,1),
                                       high=scales::muted("red"),
                                       low=scales::muted("blue"))
}
