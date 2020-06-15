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
umap <- function( X, excl )
{
    trm <- intersect( excl, colnames(Y) )
    U <- X %>% dplyr::select( -dplyr::one_of(trm) ) %>% uwot::umap() %>%
        as.data.frame() %>% dplyr::rename_all( stringr::str_replace, "V", "UMAP" ) %>%
        dplyr::mutate_all( ~(.x - min(.x))/(max(.x) - min(.x)) )
    dplyr::bind_cols(X, U)
}
