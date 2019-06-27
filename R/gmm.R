## Gaussian mixture modeling
##
## by Artem Sokolov

## Modifies a vector to have a single "activation" phase
singleAct <- function( v )
{
    n <- length(v)
    mni <- which.min(v)
    mxi <- n - which.max(rev(v)) + 1
    stopifnot( mni < mxi )
    v[1:mni] <- min(v)
    v[mxi:n] <- max(v)
    v
}

## Fits a two-component mixture model to the provided 1-D vector
## v - vector of values
## chName - name of the channel (for reporting only)
## rs - random seed to allow reproducibility
mixEM <- function( v, chName, rs=100 )
{
    set.seed(rs)
    cat( "Fitting a mixture model to", chName, "... " )
    mix <- mixtools::normalmixEM( v, mu=c(-0.8,0.8), maxit=2000 )
    mix[c("lambda","mu","sigma")]
}

## Given a model returned by mixEM(), computes posterior probabilities
##   for column vals in data frame .df
mutate_probs <- function( .df, vals, mix )
{
    ## Ensure the means are ordered
    stopifnot( mix$mu[1] < mix$mu[2] )
    s <- ensym(vals)

    ## Compute posterior probabilities of each value landing in negative and positive clusters
    .df %>% arrange(!!s) %>%
        mutate( CN = mix$lambda[1]*dnorm(!!s, mix$mu[1], mix$sigma[1]),
               CP = mix$lambda[2]*dnorm(!!s, mix$mu[2], mix$sigma[2]),
               RawProb = CP/(CP+CN), Pos = singleAct(RawProb), Neg = 1-Pos )
}

## Fit a Gaussian model to a set of channels
## X - matrix of marker expression
## qq - value between 0 and 1. For each channel, the qq^th quantile will be mapped to 0 and (1-qq)^th quantile to 1.
## ... - columns to fit a GMM to
## seed - random seed to allow for reproducibility
fitGMM <- function(X, qq, ..., seed=100)
{
    ## Argument verification
    if( qq < 0 | qq > 0.5 ) stop( "qq argument must be in [0, .5] range" )
    
    ## Isolate the marker values of interest
    MV <- X %>% gather( Marker, Values, ... ) %>% select( Marker, Values ) %>%
        group_by( Marker ) %>% summarize_at("Values",list)
    
    ## Verify that the data has been log-normalized
    if( range(MV$Values)[2] > 1000 )
        warning( "Large values detected. Please ensure the data has been log-normalized." )

    ## Compute qq^th and (1-qq)^th quantiles
    ## Adjust and filter the values accordingly
    MVQ <- MV %>% mutate( QQ = map( Values, ~list(lo=quantile(.x, qq), hi=quantile(.x, 1-qq)) ) ) %>%
        mutate_at( "Values", map2, .$QQ, ~((.x-.y$lo) / (.y$hi-.y$lo)) ) %>%
        mutate_at( "Values", map, keep, ~(.x >= 0 & .x <=1) )
    
    ## Fit a mixture of two Gaussians to each marker
    MVQ %>% mutate( GMM = map2(Values, Marker, mixEM, seed) ) %>%
        mutate_at( "GMM", map2, .$QQ, c ) %>% select( Marker, GMM )
    
    ## P <- X %>% nest( -Marker, .key=Values ) %>% inner_join( MX, by="Marker" ) %>%
    ##     mutate( PostProb = map2(Values, Fit, ~mutate_probs(.x, "Value", .y)) ) %>%
    ##     select( -Values, -Fit ) %>% unnest()
}

## Compute the probability of marker expression, using previously trained GMM models
## X - matrix of marker expression values
##probGMM <- function(X, GMMs, 
