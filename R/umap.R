## Functions related to UMAP plot generation
##
## by Artem Sokolov

## Generates a palette for a vector cell type labels
##
## Not exported
makePal <- function( v )
{
    pal <- c("Medium", "Dark", "Light") %>%
        purrr::map( ggthemes::few_pal ) %>%
        purrr::map( ~.(8) ) %>% unlist
    
    v1 <- setdiff( v, "Other (None)" )
    set_names( pal[1:length(v1)], v1 ) %>%
        c( "Other (None)" = "gray" )
}

#' Augments a matrix of marker probabilities with its UMAP projection
#' @param X data frame of probabilities or raw expression
#' @param excl character vector of column names to exclude from modeling
#' @return X, augmented with columns UMAP1, UMAP2
#' @importFrom magrittr %>%
#' @export
umap <- function( X, excl )
{
    trm <- intersect( excl, colnames(Y) )
    U <- X %>% dplyr::select( -dplyr::one_of(trm) ) %>% uwot::umap() %>%
        as.data.frame() %>% dplyr::rename_all( stringr::str_replace, "V", "UMAP" ) %>%
        dplyr::mutate_all( ~(.x - min(.x))/(max(.x) - min(.x)) )
    dplyr::bind_cols(X, U)
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
        pal <- makePal( Z$Label )
        
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


