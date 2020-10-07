##### load libraries #####

library( altR2 )
library( PearsonDS )


##### create the design dataframe #####

# factor levels: sample size 'N', number of predictors 'p', population multiple correlation
# coefficient 'rho2' and distribution 'dist'
N <- c( 20, 40, 60, 100, 150 )
p <- c( 2, 5, 10 )
rho2 <- c( 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9 )
dist <- c( 'conm', 'conh', 'bin', 'ord5', 'ord7', 'mix2', 'mix5', 'mix10' )

# design dataframe
design <- expand.grid( N, p, rho2, dist )
colnames( design ) <- c( 'N', 'p', 'rho2', 'dist' )

# remove the redundant p's in the mixed distributions, as each mixed distribution only
# consists of one factor level of the number of predictors factor
design0 <- design[ design$dist != 'mix2' &
                   design$dist != 'mix5' &
                   design$dist != 'mix10', ]
mix2 <- design[ design$dist == 'mix2' & design$p == 2, ]
mix5 <- design[ design$dist == 'mix5' & design$p == 5, ]
mix10 <- design[ design$dist == 'mix10' & design$p == 10, ]
design <- rbind( design0, mix2, mix5, mix10 )
row.names( design ) <- 1 : 630


##### data generation function #####

# input: factor level values
# output: a list containing a vector of outcome variable values 'y' and a matrix of independent
# variable values 'X'
datagen <- function( N, p, rho2, dist ) {
  
  # constants
  mean <- 0                                           # mean of continuous variables
  var <- 1                                            # variance of continuous variables
  skewm <- 0.5                                        # moderate skewness parameter
  skewh <- 1.5                                        # high skewness parameter
  kurtm <- 5                                          # moderate kurtosis parameter
  kurth <- 7                                          # high kurtosis parameter
  xbin <- c( 1, 2 )                                   # binary variable values
  pbin <- c( 2/3, 1/3 )                               # probability (prob) of binary variables
  vbin <- pbin[ 2 ] * pbin[ 1 ]                       # variance (var) of binary variables
  xord5 <- c( 1 : 5 )                                 # ordinal variable values (5 categories)
  pord5 <- c( .1, .3, .2, .3, .1 )                    # prob of ordinal variables (5 categories)
  vord5 <- sum( xord5^2 * pord5 ) - sum( xord5 * pord5 )^2                          
                                                      # var of ordinal variables (5 categories)
  xord7 <- c( 1 : 7 )                                 # ordinal variable values (7 categories)
  pord7 <- c( 1/14, 2/14, 3/14, 2/14, 3/14, 2/14, 1/14 )
                                                      # prob of ordinal variables (7 categories)
  vord7 <- sum( xord7^2 * pord7 ) - sum( xord7 * pord7 )^2
                                                      # var of ordinal variables (7 categories)
  momentsm <- c( mean, var, skewm, kurtm )            # Pearson moments (moderate parameters)
  momentsh <- c( mean, var, skewh, kurth )            # Pearson moments (high parameters)
  ones <- rep( 1, p )                                 # vector of ones (used for later equations)
  
  # the generation of the data is specified per distribution within the design:
  
  # continuous distribution with moderate levels of skewness and kurtosis
  if( dist == 'conm' ) {
    X <- matrix( 9999, N, p )
    for( i in 1 : p ) {
      X[ , i ] <- rpearson( N, moments = momentsm )
    }
    vars <- rep( var, p )
  }
  
  # continuous distribution with high levels of skewness and kurtosis
  if( dist == 'conh' ) {
    X <- matrix( 9999, N, p )
    for( i in 1 : p ) {
      X[ , i ] <- rpearson( N, moments = momentsh )
    }
    vars <- rep( var, p )
  }
  
  # binomial distribution
  if( dist == 'bin' ) {
      X <- matrix( 9999, N, p )
      for( i in 1 : p ) {
        X[ , i ] <- sample( xbin, N, prob = pbin, replace = TRUE )
      }
      vars <- rep( vbin, p )
  } 
  
  # distribution with five-class ordered categorical variables, treated as numerical
  if( dist == 'ord5' ) {
    X <- matrix( 9999, N, p )
    for( i in 1 : p ) {
      X[ , i ] <- sample( xord5, N, prob = pord5, replace = TRUE )
    }
    vars <- rep( vord5, p )
  }
  
  # distribution with seven-class ordered categorical variables, treated as numerical
  if( dist == 'ord7' ) {
    X <- matrix( 9999, N, p )
    for( i in 1 : p ) {
      X[ , i ] <- sample( xord7, N, prob = pord7, replace = TRUE )
    }
    vars <- rep( vord7, p )
  }
  
  # distribution with mixed variable types: one binary, one continuous (moderate parameters)
  if( dist == 'mix2' ) {
      x1 <- sample( xbin, N, prob = pbin, replace = TRUE )
      x2 <- rpearson( N, moments = momentsh )
      X <- cbind( x1, x2 )
    vars <- c( vbin, var )
    ones <- c( 1, 1 )
  }
  
  # distribution with mixed variable types: two binary, two continuous (one with moderate
  # parameters, one with high parameters) and one five-class ordered categorical
  if( dist == 'mix5' ) {
    x1 <- sample( xbin, N, prob = pbin, replace = TRUE )
    x2 <- sample( xbin, N, prob = pbin, replace = TRUE )
    x3 <- sample( xord5, N, prob = pord5, replace = TRUE )
    x4 <- rpearson( N, moments = momentsm )
    x5 <- rpearson( N, moments = momentsh )
    X <- cbind( x1, x2, x3, x4, x5 )
    vars <- c( vbin, vbin, vord5, var, var )
    ones <- rep( 1, 5 )
  }
  
  # distribution with mixed variable types: three binary, three continuous (one with moderate
  # parameters, two with high parameters), two five-class ordered categorical and two
  # seven-class ordered categorical
  if( dist == 'mix10' ) {
      x1 <- sample( xbin, N, prob = pbin, replace = TRUE )
      x2 <- sample( xbin, N, prob = pbin, replace = TRUE )
      x3 <- sample( xbin, N, prob = pbin, replace = TRUE )
      x4 <- sample( xord5, N, prob = pord5, replace = TRUE )
      x5 <- sample( xord5, N, prob = pord5, replace = TRUE )
      x6 <- sample( xord7, N, prob = pord7, replace = TRUE )
      x7 <- sample( xord7, N, prob = pord7, replace = TRUE )
      x8 <- rpearson( N, moments = momentsm )
      x9 <- rpearson( N, moments = momentsh )
      x10 <- rpearson( N, moments = momentsh )
      X <- cbind( x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 )
    vars <- c( rep( vbin, 3 ), vord5, vord5, vord7, vord7, rep( var, 3 ) )
    ones <- rep( 1, 10 )
  }
  
  # determining the y variable values
  covmatrix <- diag( vars )
  yvar <- ( t( ones ) %*% covmatrix %*% ones ) / rho2
  evar <- yvar - rho2 * yvar
  y <- 100 + X %*% ones + rnorm( N, sd = sqrt( evar ) )
  
  return( list( y = y, X = X ) )
  
}

##### create and fill the results dataframe #####

# create the results dataframe
reps <- 1000
resultslong <- design[ rep( seq_len( nrow( design ) ), each = reps ), ]
rownames( resultslong ) <- 1 : 630000
resultslong[ , 5 : 42 ] <- 9999
estnames <- c( 'R2', 'S', 'E', 'W', 'OP1', 'OP2', 'OP5', 'P', 'C', 'OPE', 'S+', 'E+', 'W+',
               'OP1+', 'OP2+', 'OP5+', 'P+', 'C+', 'OPE+' )
colnames( resultslong )[ 5 : 23 ] <- paste( estnames, 'bias', sep = '_' )
colnames( resultslong )[ 24 : 42 ] <- paste( estnames, 'MSE', sep = '_' )

# fill the results dataframe
set.seed( 80085 )

for( i in 1 : 630000 ) {
  N <- resultslong$N[ i ]
  p <- resultslong$p[ i ]
  rho2 <- resultslong$rho2[ i ]
  dist <- resultslong$dist[ i ]
    while( TRUE ) {
      try({
        sim <- datagen( N, p, rho2, dist )
        model <- lm( sim$y ~ sim$X )
        ests <- as.vector( altR2( model )[ -11 ] )
        bias <- ests - rho2
        MSE <- ( ests - rho2 )^2
        resultslong[ i, 5 : 42 ] <- c( bias, MSE )
      })
      if( !is( model, 'try-error' ) ) break
    }
  cat( paste( i, '\n' ) )
}




