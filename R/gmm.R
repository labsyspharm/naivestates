## Gaussian mixture modeling
##
## by Artem Sokolov

## Modifies a vector to have a single "activation" phase
##
## Not exported
singleAct <- function( v )
{
    ## Identify the "activation" phase by its positive slope
    n <- length(v)
    ph <- which( (v[2:n] - v[1:(n-1)]) > 0 )
    mni <- min(ph)
    mxi <- max(ph) + 1
    stopifnot( mni < mxi )
    v[1:mni] <- min(v)
    v[mxi:n] <- max(v)
    v
}

## Linear rescale of a vector to the [lo, hi] range
##
## Not exported
rescaleLin <- function( v, lo, hi )
{
    vadj <- (v-lo) / (hi-lo)
    vfin <- vadj[ is.finite(vadj) ]
    dplyr::case_when(
               vadj == -Inf ~ min(vfin),
               vadj == Inf ~ max(vfin),
               TRUE ~ vadj )
}

## Automatically finds distribution boundaries that contain most of the data
##
## Not exported
findBounds <- function( .v, baseline=1e-3 )
{
    density( .v ) %>% with( tibble::tibble(x, y) ) %>%
        dplyr::mutate_at( "y", ~.x / max(.x) ) %>%
        dplyr::arrange(x) %>% dplyr::filter( y > baseline ) %>%
        dplyr::filter( dplyr::row_number() %in% c(1,n()) ) %>%
        dplyr::pull(x) %>% as.list() %>% rlang::set_names(c("lo","hi"))
}

## Fits a two-component mixture model to the provided 1-D vector
## v - vector of values
## chName - name of the channel (for reporting only)
## mu_init - vector of two values, each between 0 and 1, specifying the initial placement of Gaussian means
## rs - random seed to allow reproducibility
##
## Not exported
mixEM <- function( v, chName, mu_init, rs=100 )
{
    set.seed(rs)
    cat( "Fitting a mixture model to", chName, "... " )
    mix <- mixtools::normalmixEM( v, mu=mu_init, maxit=2000 )

    ## Reorder the means to be monotonically increasing
    mdl <- mix[c("lambda","mu","sigma")]
    j <- order(mdl$mu)
    purrr::map( mdl, ~.x[j] )
}

## Given a model returned by mixEM(), computes posterior probabilities
##   for column vals in data frame .df
##
## Not exported
mutate_probs <- function( .df, vals, mix )
{
    ## Ensure the means are ordered
    stopifnot( mix$mu[1] < mix$mu[2] )
    s <- rlang::ensym(vals)

    ## Adjust the values based on quantile limits
    ## Compute posterior probabilities of each value landing in negative and positive clusters
    ## NAs arising due to CP+CN == 0 are filled, NAs in column !!s are then restored
    .df %>% dplyr::arrange(!!s) %>%
        dplyr::mutate( AdjVal = rescaleLin(!!s, mix$lo, mix$hi),
                      CN = mix$lambda[1]*dnorm(AdjVal, mix$mu[1], mix$sigma[1]),
                      CP = mix$lambda[2]*dnorm(AdjVal, mix$mu[2], mix$sigma[2]),
                      Prob = tidyr::replace_na(CP/(CP+CN), NA) ) %>%     # NaN -> NA for fill()
            tidyr::fill( Prob, .direction="downup" ) %>%
            dplyr::mutate( Prob = singleAct(Prob) ) %>%
            dplyr::mutate( Prob = ifelse(is.na(!!s), NA, Prob) )   # Restore original NAs
}

#' Fit a Gaussian model to a set of channels
#' 
#' @param X - cell-by-channel data frame of marker expression
#' @param cid - column name or index that contains cell IDs
#' @param ... - columns to fit a GMM to
#' @param bounds - [optional] Two values specifying the range of values to be used
#'             for mixture modeling. This allows for exclusion of outliers from
#'             model fitting. Calculated automatically, if not provided.
#' @param baseline - [optional] baseline for automatic boundary calculation
#' @param mu_init - [optional] vector of two values, each between 0 and 1,
#'                  specifying the initial placement of Gaussian means.
#' @param seed - [optional] random seed to allow for reproducibility
#' @return A composite data frame containing trained GMMs and their fits to data
#' @importFrom magrittr %>%
#' @export
GMMfit <- function(X, cid, ..., bounds, baseline=0.01, mu_init=c(0.2,0.8), seed=100)
{
    ## Argument verification
    mb <- missing(bounds)
    if( !mb && length(bounds) != 2 )
        stop( "bounds must be a vector of two values" )
    if( any(mu_init < 0) | any(mu_init > 1) )
        stop( "mu_init argument values must be in [0, 1] range" )

    ## Select all relevant columns from the original data frame
    X1 <- X %>% dplyr::select( {{cid}}, ... ) %>%
        tidyr::gather( Marker, Values, -1 )

    ## Verify that the data has been log-normalized
    if( range(X1$Values, na.rm=TRUE)[2] > 1000 )
        warning( "Large values detected. Please ensure the data has been log-normalized." )

    ## Isolate the finite marker values to use for modeling
    MV <- X1 %>% dplyr::filter( is.finite(Values) ) %>%
        dplyr::group_by( Marker ) %>%
        dplyr::summarize_at("Values",list)

    ## Computes distribution bounds, or uses provided values
    ## Adjust and filter the values accordingly
    fb <- ~set_names(`if`( mb, findBounds(.x,baseline), as.list(sort(bounds)) ),
                     c("lo","hi"))
    MVQ <- MV %>% dplyr::mutate( Bounds = purrr::map(Values, fb) ) %>%
        dplyr::mutate_at( "Values", map2, .$Bounds, ~((.x-.y$lo) / (.y$hi-.y$lo)) ) %>%
        dplyr::mutate_at( "Values", map, keep, ~(.x >= 0 & .x <=1) )
    
    ## Fit a mixture of two Gaussians to each marker
    G <- MVQ %>%
        dplyr::mutate( GMM = purrr::map2(Values, Marker, mixEM, mu_init, seed) ) %>%
        dplyr::mutate_at( "GMM", purrr::map2, .$Bounds, c ) %>%
        dplyr::select( Marker, GMM )
    
    ## Compute posterior probabilities on the original data
    fmp <- ~mutate_probs(.x, "Value", .y)
    R <- X1 %>% dplyr::rename( Value = Values ) %>%
        tidyr::nest( Values=c({{cid}}, Value) ) %>%
        dplyr::inner_join(G, by="Marker") %>%
        dplyr::mutate( Values = purrr::map2(Values, GMM, fmp) )
    R
}

#' Reshapes a GMMfit() dataframe to the original cell-by-marker format
#' 
#' @param .df - data frame produced by GMMfit()
#' @return A data frame in the original cell-by-marker format
#' @importFrom magrittr %>%
#' @export
GMMreshape <- function(.df)
{
    .df %>% dplyr::select( -GMM ) %>% tidyr::unnest(Values) %>%
        dplyr::select( -Value, -AdjVal, -CN, -CP ) %>%
            tidyr::spread( Marker, Prob )
}

#' Returns GMM models to tidy format
#'
#' @param .df - data frame produced by GMMfit()
#' @return A data frame that maps each marker to the following values in
#'     the log10-transformed space:
#' \describe{
#'   \item{lo/hi}{The low and high cut-off values for outlier exclusion.}
#'   \item{lambda1/2}{Mixture coefficients for the GMM. Sum up to 1.}
#'   \item{mu1/2}{Centroids of the two Gaussians in the mixture model.}
#'   \item{sigma1/2}{The corresponding standard deviation of each Gaussian.}
#' }
#' @importFrom magrittr %>%
#' @export
GMMmodels <- function(.df)
{
    .df %>%
    dplyr::select( Marker, GMM ) %>%
    dplyr::mutate_at( "GMM", purrr::map,
                     ~list(lo = .x$lo, hi = .x$hi,
                           lambda1 = .x$lambda[1],
                           lambda2 = .x$lambda[2],
                           mu1 = .x$mu[1] * (.x$hi - .x$lo) + .x$lo,
                           mu2 = .x$mu[2] * (.x$hi - .x$lo) + .x$lo,
                           sigma1 = .x$sigma[1] * (.x$hi - .x$lo),
                           sigma2 = .x$sigma[2] * (.x$hi - .x$lo))
                     ) %>%
    dplyr::mutate_at( "GMM", purrr::map,
                     ~tidyr::unnest(tibble::enframe(.x),value) ) %>%
    tidyr::unnest(GMM) %>% tidyr::pivot_wider()
}
