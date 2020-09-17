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
