score.crossing <- function(weight.1, weight.2, N.cross){
        
    # determine the index of the split branch in weight.1
    parent <- which( !(rownames(weight.1) %in% rownames(weight.2) ) )
    
    # determine the indices of the children in weight.2
    children <- which( !(rownames(weight.2) %in% rownames(weight.1) ) )

    # check the provided matrices
    if (length(parent) != 1) {
    stop("Wrong weight matrices provided. The rownames do not match...")}
    check <- any(weight.1[-parent, ] != weight.2[-children, ])
    if (check) {stop("The matrices provided do not correspond to two trees
        diferring in one expanded branch...")}
        
    # compute the score before splitting
    size.hierarchical.1 <- apply(weight.1, 1, sum)
    size.flat <- apply(weight.1, 2, sum)
    jaccard.1 <- sum(weight.1[parent,] / (size.hierarchical.1[parent] +
        size.flat - weight.1[parent, ]))
    common.1 <- sum(weight.1[parent, ]^2)
    score.1 <- jaccard.1 * common.1
    
    # compute the score after splitting. Note: the crossing provided must
    # correspond to the specific layout of the bigraph
    size.hierarchical.2 <- apply(weight.2, 1, sum)
    jaccard.2 <- sum(weight.2[children, ] /
        sweep(size.hierarchical.2[children] -
        weight.2[children, ], 2, size.flat, "+"))
    common.2 <- sum(weight.2[children, ]^2)
    score.2 <- jaccard.2 * common.2 - N.cross
    
    return(list(sc1 = score.1, sc2 = score.2))
}