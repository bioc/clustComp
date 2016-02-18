score.it <- function(weight.1, weight.2){

    # determine the index of the split branch in weight.1
    parent <- which( !(rownames(weight.1) %in% rownames(weight.2) ) )

    # determine the indices of the children in weight.2
    children <- which( !(rownames(weight.2) %in% rownames(weight.1) ) )
    
    # check the provided matrices
    if (length(parent) != 1) {
        stop("Wrong weight matrices provided: the rownames do not match")}
    check <- any(weight.1[-parent, ] != weight.2[-children, ])

    if (check) {
        stop("The matrices provided do not correspond to two trees diferring in
            one expanded branch...")
    } # end IF

    # compute the score before splitting
    size.hierarchical.1 <- apply(weight.1, 1, sum)
    size.flat <- apply(weight.1, 2, sum)
    p.parent.ij.A <- weight.1[parent, ] / size.hierarchical.1[parent]

    # vector adding up to 1
    p.parent.ij.B <- weight.1[parent, ] / size.flat
    index <- which(weight.1[parent,] != 0)
    N <- sum(weight.1)
    score.1 <- -sum(log2( p.parent.ij.A[index] * p.parent.ij.B[index] ) *
        (weight.1[parent, index] + 1)) / (2 * N)

    # compute the score after splitting
    size.hierarchical.2 <- apply(weight.2, 1, sum)
    p.children.ij.A <- weight.2[children, ] / size.hierarchical.2[children]
    p.children.ij.B <- sweep(weight.2[children, ], 2, size.flat, "/")
    index <- which(weight.2[children, ] != 0, arr.ind = TRUE)
    score.2 <- -sum(log2( p.children.ij.A[index] * p.children.ij.B[index]) *
        (weight.2[children, ][index] + 1)) / (2 * N)

    return(list(sc1 = score.1, sc2 = score.2))
}
