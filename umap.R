## Generates a palette for a vector cell type labels
makePal <- function( v )
{
    pal <- c("Medium", "Dark", "Light") %>%
        map( ggthemes::few_pal ) %>%
        map( ~.(8) ) %>% unlist
    
    v1 <- setdiff( v, "Other (None)" )
    set_names( pal[1:length(v1)], v1 ) %>%
        c( "Other (None)" = "gray" )
}

## A UMAP plot summarizing cell state calls
plotSummary <- function( Y )
{
    cat( "Computing UMAP projection...\n" )
    
    ## Compute a UMAP projection of the data
    trm <- intersect( c(opt$id, "State", "Anchor"), colnames(Y) )
    U <- Y %>% select( -one_of(trm) ) %>% uwot::umap() %>%
        as.data.frame() %>% rename_all( str_replace, "V", "UMAP" ) %>%
        mutate_all( ~(.x - min(.x))/(max(.x) - min(.x)) )

    ## Augment the original data
    Z <- bind_cols(Y, U)

    ## Color by cell state calls (if any)
    if( "State" %in% colnames(Z) ) {
        Z <- Z %>% mutate( Label = str_c(State, " (", Anchor, ")") )
        pal <- makePal( Z$Label )

        ## Plot UMAP projection, coloring by cell state calls
        ggplot( Z, aes(UMAP1, UMAP2, color=Label) ) +
            geom_point() + theme_bw() +
            scale_color_manual( values=pal, name="Cell State Call\n(Dominant Marker)" )
    } else
        ggplot( Z, aes(UMAP1, UMAP2) ) + theme_bw() + geom_point(color="gray")
}

