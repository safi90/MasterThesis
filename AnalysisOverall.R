# create a vector with full estimator names, which will be used in tables and such
estnames <- c( 'R-squared', 'Smith', 'Ezekiel', 'Wherry', 'Olkin-Pratt, K = 1',
               'Olkin-Pratt, K = 2', 'Olkin-Pratt, K = 5', 'Pratt', 'Claudy',
               'Exact Olkin-Pratt', 'Positive-part Smith', 'Positive-part Ezekiel',
               'Positive-part Wherry', 'Positive-part Olkin-Pratt, K = 1',
               'Positive-part Olkin-Pratt, K = 2', 'Positive-part Olkin-Pratt, K = 5',
               'Positive-part Pratt', 'Positive-part Claudy',
               'Positive-part exact Olkin-Pratt' )


##### bias analysis #####

# create a short results dataframe (condense repetitions to a single research design cell)
resultsshort <- design
resultsshort[ , 5 : 42 ] <- 9999
colnames( resultsshort )[ 5 : 23 ] <- paste( estnames, 'bias', sep = '_' )
colnames( resultsshort )[ 24 : 42 ] <- paste( estnames, 'MSE', sep = '_' )

for ( i in 1 : nrow( resultsshort ) ) {
  onerow <- resultslong[ resultslong$N == design$N[ i ] &
                           resultslong$p == design$p[ i ] &
                           resultslong$rho2 == design$rho2[ i ] &
                           resultslong$dist == design$dist[ i ], ]
  resultsshort[ i, 5 : 42 ] <- colMeans( onerow[ , 5 : 42 ] )
}

# get the p-values; testing the biases against zero
pvals <- resultsshort
pvals[ , 5 : 42 ] <- 9999
for( i in 1 : nrow( pvals ) ) {
  onerow <- resultslong[ resultslong$N == design$N[ i ] &
                           resultslong$p == design$p[ i ] &
                           resultslong$rho2 == design$rho2[ i ] &
                           resultslong$dist == design$dist[ i ], ]
  for( j in 5 : ncol( pvals ) )
    pvals[ i, j ] <- t.test( onerow[ , j ] )$p.value
}

# create a data frame displaying unbiased research design cells per estimator in absolute
# values and percentages
nobias <- data.frame( cbind( estnames, rep( 9999, length( estnames ) ) ) )
colnames( nobias ) <- c( 'Estimator', 'Number of cells with bias' )
nobias$`Number of cells with bias` <- as.numeric( nobias$`Number of cells with bias` )
for( i in 1 : nrow( nobias ) ) {
  ps <- pvals[ , i + 4 ]
  nobias[ i, 2 ] <- length( which( ps < 0.005 ) )
}
nobias <- nobias[ order( nobias$`Number of cells with bias`, decreasing = TRUE ), ]
for( i in 1 : nrow( nobias ) ) {
  nobias[ i, 3 ] <- round( ( nobias$`Number of cells with bias`[ i ] / sum(
    nobias$`Number of cells with bias`) * 100 ), 2 )
}
colnames( nobias )[ 3 ] <- 'Percentage'


# histogram of the overall bias results
pal <- colorRampPalette( colors = c( 'blue', 'lightblue' ) ) ( 19 )
par( mar = c( 5, 11, 4, 2 ) )
bp <- barplot( nobias$`Number of cells with bias`, names.arg = nobias$Estimator, las = 1,
               horiz = T, cex.names = 0.8, xlim = c( 0, 1500 ), col = pal )
text( y = bp, x = nobias$`Number of cells with bias`, pos = 4,
      labels = paste( nobias$`Number of cells with bias`, '(', nobias$Percentage, '%', ')' ),
      cex = 0.8 )


##### MSE dominance analysis #####

# generate a matrix with the MSE dominance relationships between the estimators
domMSE <- data.frame( matrix( 9999, length( estnames ), length( estnames ) ) )
colnames( domMSE ) <- rownames ( domMSE ) <- estnames
resultsMSE <- resultsshort[ , -c( 5 : 23 ) ]
for( i in 1 : nrow( domMSE ) ) {
  for( j in 1 : ncol( domMSE ) ) {
    if( i != j ) {
      demand <- rep( '+', nrow( design ) )
      check <- rep( 9999, nrow( design ) )
      for( row in 1 : nrow( design ) ) {
        if( resultsMSE[ row, i + 4 ] <= resultsMSE[ row, j + 4 ] ) {
          check[ row ] <- '+'
        } else {
          check[ row ] <- '-'
        }
      }
      if( all.equal( check, demand ) == TRUE ) {
        domMSE[ i, j ] <- '+'
      } else {
        domMSE[ i, j ] <- '-'
      }
    } else {
      domMSE[ i, j ] <- 'X'
    }
  }
}


##### average MSE analysis #####

# get the average MSEs
avMSE <- expand.grid( N, p )
avMSE[ , 3 : 21 ] <- 9999
colnames( avMSE ) <- c( 'N', 'p', estnames )
for( i in 1 : nrow( avMSE ) ) {
  onerow <- resultslong[ resultslong$N == avMSE$N[ i ] &
                           resultslong$p == avMSE$p[ i ], 24 : 42 ]
  avMSE[ i, 3 : 21 ] <- colMeans( onerow )
}
avMSErho2 <- expand.grid( N, p, rho2 )
avMSErho2[ , 4 : 22 ] <- 9999
colnames( avMSErho2 ) <- c( 'N', 'p', 'rho-squared', estnames )
for( i in 1 : nrow( avMSErho2 ) ) {
  onerow <- resultslong[ resultslong$N == avMSErho2$N[ i ] &
                           resultslong$p == avMSErho2$p[ i ] &
                           resultslong$rho2 == avMSErho2$`rho-squared`, 24 : 42 ]
  avMSErho2[ i, 4 : 22 ] <- colMeans( onerow )
}
avMSErho2[ , 4 : 22 ] <- round( avMSErho2[ , 4 : 22 ], 4 )

# determine which estimator has the lowest MSE in each condition (combination of sample size
# and number of predictors factor level)
lowestAMSE <- avMSE[ , 1 : 4 ]
colnames( lowestAMSE )[ c( 3, 4 ) ] <- c( 'best estimator', 'average MSE' )
for( i in 1 : nrow( lowestAMSE ) ) {
  lowestAMSE[ i, 3 ] <- colnames( avMSE )[ which.min( avMSE[ i, ] ) ]
  lowestAMSE[ i, 4 ] <- round( avMSE[ i, ][ which.min( avMSE[ i, ] ) ], 4 )
}

# expand the above with the addition of rho2 values to the factor level combinations
lowestAMSErho2 <- avMSErho2[ , 1 : 4 ]
colnames( lowestAMSErho2 )[ 4 ] <- 'best estimator'
for( i in 1 : nrow( lowestAMSErho2 ) ) {
  lowestAMSErho2[ i, 4 ] <- colnames( avMSErho2 )[ -3 ][ which.min( avMSErho2[ i, -3 ] ) ]
}

# the minimum average MSE per estimator overall
avMSEall <- rep( 9999, length( estnames ) )
for( i in 1 : length( avMSEall ) ) {
  onerho <- avMSE[ , i + 2 ]
  avMSEall[ i ] <- round( mean( onerho ), 4 )
}


##### maximum MSE analysis #####

# get the maximum MSEs
maxMSE <- avMSE
for( i in 1 : nrow( maxMSE ) ) {
  for( j in 3 : ncol( maxMSE ) ) {
    onecell <- resultsshort[ resultsshort$N == maxMSE$N[ i ] &
                               resultsshort$p == maxMSE$p[ i ], ]
    one.est <- onecell[ , j + 21 ]
    maxMSE[ i, j ] <- one.est[ which( one.est == max( one.est ) ) ]
  }
}

# identify best maximum MSE estimator per combination of sample size and number of predictors
bestworst <- maxMSE[ , 1 : 4 ]
colnames( bestworst )[ c( 3, 4 ) ] <- c( 'best estimator', 'maximum MSE' )
for( i in 1 : nrow( bestworst ) ) {
  oner <- maxMSE[ i, 3 : 21 ]
  best <- which.min( oner )
  bestworst[ i, 3 ] <- estnamesfull[ best ]
  bestworst[ i, 4 ] <- round( oner[ best ], 4 )
}

# maximum MSE of each estimator overall
maxMSEall <- rep( 9999, length( estnames ) )
for( i in 1 : length( maxMSEall ) ) {
  onerho <- maxMSE[ , i + 2 ]
  maxMSEall[ i ] <- round( onerho[ which( maxMSE[ , i + 2 ] == max( maxMSE[ , i + 2 ] ) ) ], 4 )
}

# overall average and max MSEs of the estimators
overall <- cbind( avMSEall, maxMSEall )
rownames( overall ) <- estnamesfull
colnames( overall ) <- c( 'estimator', 'average MSE', 'maximum MSE' )

