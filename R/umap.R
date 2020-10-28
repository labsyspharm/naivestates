## Functions related to UMAP plot generation
##
## by Artem Sokolov

#' Augments a matrix of marker probabilities with its UMAP projection
#' @param X data frame of probabilities or raw expression
#' @param excl character vector of column names to exclude from modeling
#' @return X, augmented with columns UMAP1, UMAP2
#' @importFrom magrittr %>%
#' @export
umap <- function( X, excl )
{
    trm <- intersect( excl, colnames(X) )
    U <- na.omit(X) %>% dplyr::select( -dplyr::one_of(trm) ) %>% uwot::umap() %>%
        as.data.frame() %>% dplyr::rename_all( stringr::str_replace, "V", "UMAP" ) %>%
        dplyr::mutate_all( ~(.x - min(.x))/(max(.x) - min(.x)) )
    dplyr::bind_cols(na.omit(X), U)
}


