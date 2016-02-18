barycentre <- function(edge.weight, coordinates=NULL){
        
    if (is.vector(edge.weight) == FALSE) {
        stop("Incorrect dimension of the edge weights...")
    }
    if (length(coordinates)==0) {
        coordinates <- c(length(edge.weight):1)
    }
    if (length(coordinates) != length(edge.weight)) {
        stop("The lengths of the edge.weight and coordinates don't match...")
    }
    if (any(edge.weight<0)) {
        stop("All weights must be non negative...")
    }
    position <- sum(coordinates*edge.weight)/sum(edge.weight)
    return(position)
}

