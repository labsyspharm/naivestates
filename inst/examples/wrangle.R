## Wrangling of all relevant files
##
## by Artem Sokolov

library( tidyverse )

## Predefined marker -> cell type mapping
markerMap <- function()
{
    tribble(
        ~Marker, ~Class,
        "CATENIN", "Tumor",
        "CD3", "Immune",
        "CD45RO", "Immune",
        "CD45", "Immune",
        "CD4", "Immune",
        "CD8", "Immune",
        "FOXP3", "Immune",
        "MITF", "Tumor",
        "PD1", "Immune",
        "S100", "Tumor",
        "SMA", "Stroma",
        "VIMENTIN", "Stroma",
        "ECAD", "Epithelial",
        "KERATIN", "Epithelial",
        "VIMENTIN", "Mesenchymal",
        "CD8", "T-ctx",
        "CD4", "T-hlp",
        "CD4", "T-reg",
        "FOXP3", "T-reg"
    )
}

## Classification task definitions
taskMap <- function()
{
    list( "Task1" = c("Tumor", "Immune", "Stroma"),
         "Task2" = c("Epithelial", "Mesenchymal"),
         "Task3" = c("T-ctx", "T-hlp", "T-reg") )
}

main <- function()
{
    ## Load the segmented cells-by-marker matrix and the marker map
    X <- read_csv( "example1_data.csv.gz", col_types=cols(CellId=col_integer()) )

    ## Compose the Channel -> Class map
    CM <- colnames(X) %>% keep( ~grepl("Cell_",.) ) %>%
        tibble( Channel=., Marker = str_sub(., 13)) %>%
        inner_join( markerMap(), by="Marker" ) %>%
        select( Channel, Class ) %>% write_csv( "example1_chnlmap.csv" )

    ## Write out the task definitions
    taskMap() %>% yaml::write_yaml( "example1_taskmap.yaml" )
}
