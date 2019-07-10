## Gaussian mixture modeling
##
## by Artem Sokolov

## Modifies a vector to have a single "activation" phase
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

## Fits a two-component mixture model to the provided 1-D vector
## v - vector of values
## chName - name of the channel (for reporting only)
## mu_init - vector of two values, each between 0 and 1, specifying the initial placement of Gaussian means
## rs - random seed to allow reproducibility
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
mutate_probs <- function( .df, vals, mix )
{
    ## Ensure the means are ordered
    stopifnot( mix$mu[1] < mix$mu[2] )
    s <- ensym(vals)

    ## Adjust the values based on quantile limits
    ## Compute posterior probabilities of each value landing in negative and positive clusters
    .df %>% arrange(!!s) %>%
        mutate( AdjVal = (!!s - mix$lo)/(mix$hi - mix$lo),
               CN = mix$lambda[1]*dnorm(AdjVal, mix$mu[1], mix$sigma[1]),
               CP = mix$lambda[2]*dnorm(AdjVal, mix$mu[2], mix$sigma[2]),
               Prob = singleAct(CP/(CP+CN)) )
}

## Fit a Gaussian model to a set of channels
## X - matrix of marker expression
## cid - column that contains cell IDs
## ... - columns to fit a GMM to
## qq - value between 0 and 0.5. For each channel, the qq^th quantile will be mapped to 0 and (1-qq)^th quantile to 1.
## mu_init - vector of two values, each between 0 and 1, specifying the initial placement of Gaussian means
## seed - random seed to allow for reproducibility
GMMfit <- function(X, cid, ..., qq=0.001, mu_init=c(0.2,0.8), seed=100)
{
    ## Argument verification
    if( qq < 0 | qq > 0.5 ) stop( "qq argument must be in [0, .5] range" )
    if( any(mu_init < 0) | any(mu_init > 1) )
        stop( "mu_init argument values must be in [0, 1] range" )
    
    ## Isolate the marker values of interest
    MV <- X %>% tidyr::gather( Marker, Values, ... ) %>%
        dplyr::select( Marker, Values ) %>%
        dplyr::group_by( Marker ) %>% dplyr::summarize_at("Values",list)
    
    ## Verify that the data has been log-normalized
    if( range(MV$Values)[2] > 1000 )
        warning( "Large values detected. Please ensure the data has been log-normalized." )

    ## Compute qq^th and (1-qq)^th quantiles
    ## Adjust and filter the values accordingly
    fq <- ~list(lo=quantile(.x, qq), hi=quantile(.x, 1-qq))
    MVQ <- MV %>% dplyr::mutate( QQ = purrr::map(Values, fq) ) %>%
        dplyr::mutate_at( "Values", map2, .$QQ, ~((.x-.y$lo) / (.y$hi-.y$lo)) ) %>%
        dplyr::mutate_at( "Values", map, keep, ~(.x >= 0 & .x <=1) )
    
    ## Fit a mixture of two Gaussians to each marker
    G <- MVQ %>%
        dplyr::mutate( GMM = purrr::map2(Values, Marker, mixEM, mu_init, seed) ) %>%
        dplyr::mutate_at( "GMM", purrr::map2, .$QQ, c ) %>%
        dplyr::select( Marker, GMM )
    
    ## Compute posterior probabilities on the original data
    fmp <- ~mutate_probs(.x, "Value", .y)
    X %>% dplyr::select( !!rlang::enquo(cid), G$Marker ) %>%
        tidyr::gather( Marker, Value, -1 ) %>%
            tidyr::nest( -Marker, .key="Values" ) %>%
            dplyr::inner_join(G) %>%
            dplyr::mutate( Values = purrr::map2(Values, GMM, fmp) )
}

## Reshapes a GMMfit() dataframe to the original cell-by-marker
## DF - data frame produced by GMMfit()
GMMreshape <- function(DF)
{
    DF %>% dplyr::select( -GMM ) %>% tidyr::unnest() %>%
        dplyr::select( -Value, -AdjVal, -CN, -CP ) %>%
            tidyr::spread( Marker, Prob )
}
