cosineClusterDist <- function(X, w1, w2){

# calcultes cosine cluster distance in the format required by slingshot 1.0.0

library( matrixStats)

    mu1 <- colWeightedMeans(X, w = w1)
    mu2 <- colWeightedMeans(X, w = w2)
    cosDist = sqrt( 2*(1 - sum(mu1*mu2)/sqrt(sum(mu1^2)*sum(mu2^2))))
    return( cosDist)
}

