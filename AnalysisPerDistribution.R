##### bias analysis #####

# create a data frame depicting number of biased research cells for each estimator in each
# distribution
biasperdist <- data.frame( matrix( 9999, 19, 9 ) )
biasperdist[ , 1 ] <- estnames
colnames( biasperdist ) <- c( 'estimator', dist )
for( i in 1 : nrow( biasperdist ) ) {
  for( j in 1 : length( dist ) ) {
    pvs <- pvals[ pvals$dist == dist[ j ], ]
    biasperdist[ i, j + 1 ] <- sum( pvs[ , i + 4 ] < 0.005 )
  }
}


##### MSE dominance analysis #####

# create a list in which to establish the MSE dominance relationships per distribution
listMD <- rep( list( 9999 ), length( dist ) )
for( iter in 1 : length( dist ) ) {
  distr <- dist[ iter ]
  dMSE <- data.frame( matrix( 9999, 19, 19 ) )
  colnames( dMSE ) <- rownames ( dMSE ) <- estnames
  resMSE <- resultsshort[ resultsshort$dist == distr, -c( 5 : 23 ) ]
  for( i in 1 : nrow( dMSE ) ) {
    for( j in 1 : ncol( dMSE ) ) {
      if( i != j ) {
        demand <- rep( '+', nrow( resMSE ) )
        check <- rep( 9999, nrow( resMSE ) )
        for( row in 1 : nrow( resMSE ) ) {
          if( resMSE[ row, i + 4 ] <= resMSE[ row, j + 4 ] ) {
            check[ row ] <- '+'
          } else {
            check[ row ] <- '-'
          }
        }
        if( all.equal( check, demand ) == TRUE ) {
          dMSE[ i, j ] <- '+'
        } else {
          dMSE[ i, j ] <- '-'
        }
      } else {
        dMSE[ i, j ] <- 'X'
      }
    }
  }
  listMD[[ iter ]] <- dMSE
  names( listMD )[ iter ] <- dist[ iter ]
}


##### average MSE analysis #####

# create a list with the average MSEs per cell per distribution
listAMSE <- rep( list( 9999 ), length( dist ) )
for( i in 1 : 8 ) {
  distr <- dist[ i ]
  AMSE <- expand.grid( N, p )
  AMSE[ , 3 : 21 ] <- 9999
  colnames( AMSE ) <- c( 'N', 'p', estnames )
  for( j in 1 : nrow( AMSE ) ) {
    onerow <- resultslong[ resultslong$N == AMSE$N[ j ] &
                             resultslong$p == AMSE$p[ j ] &
                             resultslong$dist == distr, 24 : 42 ]
    AMSE[ j, 3 : 21 ] <- colMeans( onerow )
  }
  listAMSE[[ i ]] <- na.omit( AMSE )
  names( listAMSE )[ i ] <- dist[ i ]
}
names( listAMSE ) <- dist

# establish which estimators have the lowest MSE per combination of sample size and number
# of predictions per distribution
listLAM <- rep( list( 9999 ), 8 )
for( i in 1 : 8 ) {
  lowestAM <- listAMSE[[ i ]][ , 1 : 4 ]
  colnames( lowestAM )[ c( 3, 4 ) ] <- c( 'best estimator', 'average MSE' )
  for( j in 1 : nrow( lowestAM ) ) {
    lowestAM[ j, 3 ] <- colnames( listAMSE[[ 1 ]] )[ which.min( listAMSE[[ i ]][ j, ]) ]
    lowestAM[ j, 4 ] <- listAMSE[[ i ]][ j, ][ which.min( listAMSE[[ i ]][ j, ]) ]
  }
  listLAM[[ i ]] <- na.omit( lowestAM )
  names( listLAM )[ i ] <- dist[ i ]
}

# summary of the average MSEs per distribution
avMSEdists <- data.frame( matrix( 9999, length( dist ), length( estnames ) ) )
for( i in 1 : length( dist ) ) {
  av <- colMeans( listAMSE[[ i ]][ , 3 : 21 ] )
  avMSEdists [ i, ] <- av
}
colnames( avMSEdists ) <- estnames
rownames( avMSEdists ) <- dist
bestperdist <- data.frame( matrix( 9999, length( dist ), 3 ) )
bestperdist[ , 1 ] <- dist
colnames( bestperdist ) <- c( 'distribution', 'estimator', 'average MSE' )
for( i in 1 : length( dist ) ) {
  bestperdist[ i, 2 ] <- estnames[ which( avMSEdists[ i, ] == min( avMSEdists[ i, ] ) ) ]
  bestperdist[ i, 3 ] <- round( 
    avMSEdists[ i, ][ which( avMSEdists[ i, ] == min( avMSEdists[ i, ] ) ) ], 4 )
}


##### maximum MSE analysis #####

# identify maximum estimator values per condition per distribution
listMM <- rep( list( 9999 ), length( dist ) )
for( iter in 1 : length( dist ) ) {
  distr <- dist[ iter ]
  if( distr == 'mix2' ) {
    MM <- expand.grid( N, 2 )
  } else if( distr == 'mix5' ) {
    MM <- expand.grid( N, 5 )
  } else if( distr == 'mix10' ) {
    MM <- expand.grid( N, 10 )
  } else {
    MM <- expand.grid( N, p )
  }
  MM[ , 3 : 21 ] <- 9999
  colnames( MM ) <- c( 'N', 'p', estnames )
  for( i in 1 : nrow( MM ) ) {
    for( j in 3 : ncol( MM ) ) {
      onecell <- resultsshort[ resultsshort$N == MM$N[ i ] &
                                 resultsshort$p == MM$p[ i ] &
                                 resultsshort$dist == distr, ]
      one.est <- onecell[ , j + 21 ]
      max <- which( one.est == max( one.est ) )
      MM[ i, j ] <- one.est[ max ]
    }
  }
  listMM[[ iter ]] <- MM
  names( listMM )[[ iter ]] <- dist[ iter ]
}

# summary of maximum MSE per distribution
maxMSEdists <- data.frame( matrix( 9999, length( dist ), length( estnames ) ) )
colnames( maxMSEdists ) <- estnames
rownames( maxMSEdists ) <- dist
for( i in 1 : length( dist ) ) {
  max <- rep( 9999, length( estnames ) )
  for( j in 1 : length( estnames ) ) {
    max[ j ] <- listMM[[ i ]][ , j + 2 ][ which( listMM[[ i ]][ , j + 2 ] == max(
      listMM[[ i ]][ , j + 2 ] ) ) ]
  }
  maxMSEdists[ i, ] <- max
}
bpdMM <- data.frame( matrix( 9999, 8, 3 ) )
bpdMM[ , 1 ] <- dist
colnames( bpdMM ) <- c( 'distribution', 'estimator', 'maximum MSE' )
for( i in 1 : length( dist ) ) {
  bpdMM[ i, 2 ] <- estnamesfull[ which( maxMSEdists[ i, ] == min( maxMSEdists[ i, ] ) ) ]
  bpdMM[ i, 3 ] <- round(
    maxMSEdists[ i, ][ which( maxMSEdists[ i, ] == min( maxMSEdists[ i, ] ) ) ], 4 )
}