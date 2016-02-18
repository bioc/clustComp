dyn.cross <- function(weights){

    crossings <- 0
    weights <- as.matrix(weights)
    m <- nrow(weights); n <- ncol(weights)

    if ((m != 1) & (n != 1)) {
        weights <- weights[, n:1]
        c <- matrix(0, m, n)
        for (i in 2:m){
            for (j in 2:n){
                c[i,j] <- c[i, j - 1] + c[i - 1, j] - c[i - 1, j - 1] +
                    weights[i - 1, j - 1]
                crossings <- crossings + c[i, j] * weights[i, j]
            } # end FOR j
        } # end FOR i
    } # end IF

    return(crossings)
}