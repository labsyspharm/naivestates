## Marker - cell type associations
##
## by Artem Sokolov

#' Finds markers among a vector of names
#'
#' @param v - vector of names to examine
#' @param mrk - vector of markers to search for
#' @param errOnNotFound - what to do when a requested marker is not found:
#'             TRUE - throw an error; FALSE (default) - exclude from results
#' @param errOnNonUnique - what to do when a mapping is not unique:
#'             TRUE - throw an error; FALSE (default) - return the first match
#' @return A named vector of matches between v and mrk
#' @importFrom magrittr %>%
#' @export
findMarkers <- function(v, mrk, errOnNotFound=FALSE, errOnNonUnique=FALSE)
{
    ## Search for markers among the candidate names
    mrki <- mrk %>% rlang::set_names() %>%
        purrr::map_chr( stringr::str_c, "$" ) %>%
        purrr::map( grep, v )

    ## Handle missing matches
    if( errOnNotFound )
    purrr::iwalk( mrki, ~if(length(.x) == 0)
                             stop("Marker ", .y, " is not found") )
    mrki <- purrr::keep( mrki, ~length(.x) > 0 )
    
    ## Handle multiple hits
    if( errOnNonUnique )
    purrr::iwalk( mrki, ~if(length(.x) > 1)
                             stop("Marker ", .y, " maps to multiple columns") )
    mrki <- purrr::map_if( mrki, ~length(.x) > 1, ~.x[1] )

    ## Index the original vector
    purrr::map_chr( mrki, ~v[.x] )
}

#' Automatically finds markers by removing blacklisted items from a vector of names
#'
#' @param v vector of names
#' @return A subset of v that is not blacklisted
#' @export
autoMarkers <- function( v )
{
    omit <- c("AF488", "AF555", "AF647", "A488", "A555", "A647", "DNA",
              "X_position", "Y_position",
              "X_centroid", "Y_centroid", "column_centroid", "row_centroid", 
              "Area", "MajorAxisLength", "MinorAxisLength", "Eccentricity", 
              "Solidity", "Extent", "Orientation")

    j <- purrr::map( omit, grep, v ) %>% unlist
    v[setdiff( 1:length(v), j )]
}
