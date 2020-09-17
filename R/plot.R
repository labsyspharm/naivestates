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
    X1 <- FT %>% dplyr::filter( Marker == marker ) %>% purrr::pluck("Values",1) %>%
        na.omit() %>% dplyr::filter( AdjVal >= 0, AdjVal <= 1 ) %>% dplyr::arrange(Value)
    
    gg1 <- ggplot2::ggplot( X1, aes(x=AdjVal) ) + ggplot2::theme_bw() +
        bold_theme() +
        ggplot2::ylab("Density") + ggplot2::xlab("Normalized Expression") +
        ggplot2::geom_density(lwd=1.1) +
        ggplot2::facet_wrap( ~str_c(marker) )
    if( plotFit )
        gg1 <- gg1 +
            ggplot2::geom_line( aes(y=CP), color="red", lwd=1.1 ) +
            ggplot2::geom_line( aes(y=CN), color="blue", lwd=1.1 )
    
    if( !plotMap ) return(gg1)
    
    gg2 <- ggplot2::ggplot( X1, aes(x=AdjVal, y=Prob) ) +
        ggplot2::geom_line(lwd=1.1) + ggplot2::theme_bw() +
        bold_theme() +
        ggplot2::ylab( "P( Expressed )" ) + ggplot2::xlab("Normalized Expression") +
        ggplot2::scale_y_continuous( breaks = seq(0,1,by=0.2), limits=c(0,1) )
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
    vs <- 1:nrow(FT[[1]])
    if( length(vs) > npt ) vs <- sample(vs, npt)

    X <- FT %>% dplyr::mutate_at( "Values", map, slice, vs ) %>%
        dplyr::select( Marker, Values ) %>% tidyver::unnest( Values ) %>%
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
    if( all(v) && length(unique(U$Anchor)) <= 24 )
    {
        U <- dplyr::mutate( U, Label = stringr::str_c(State, " (", Anchor, ")") )
        pal <- makePal( U$Label )
        
        ## Plot UMAP projection, coloring by cell state calls
        ggplot2::ggplot( U, ggplot2::aes(UMAP1, UMAP2, color=Label) ) +
            ggplot2::geom_point() + ggplot2::theme_bw() +
                ggplot2::scale_color_manual( values=pal,
                                            name="Cell State Call\n(Dominant Marker)" )
    }
    else
        ggplot2::ggplot( U, aes(UMAP1, UMAP2) ) +
            ggplot2::theme_bw() + ggplot2::geom_point(color="gray")
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
                                       high=scales::muted("red"),
                                       low=scales::muted("blue"))
}
