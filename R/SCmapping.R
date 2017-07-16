SCmapping <- function(clustering1, clustering2, plotting = TRUE, h.min = 0.1,
    line.wd = 3, point.sz = 3, offset = 0.1, evenly = TRUE, horiz = FALSE,
    max.iter = 24, node.col = NULL, edge.col = NULL, ...){

    # initial setup and checkings
    if (length(clustering1) != length(clustering2))
        {stop("The lengths of both clusterings do not match...")}
    m <- length(unique(clustering1)); n <- length(unique(clustering2))
    weights <- table(clustering1, clustering2)
    zero.one <- matrix(0, m, n)
    colnames(zero.one) <- colnames(weights)
    rownames(zero.one) <- rownames(weights)
    myMax1 <- apply(weights, 1, max)
    myMax2 <- apply(weights, 2, max)

    ## greedy by rows
    for (a.i in 1:m) {
        merge <- which(weights[a.i,] == myMax1[a.i])
        zero.one[a.i, merge] <- 1
    } ## use FOR loop to find ALL clusters with maxOverlap

    ## greedy by cols
    for (b.j in 1:n){
        merge <- which(weights[, b.j] == myMax2[b.j])
        zero.one[merge,b.j] <- 1
    }

    ## fill in with ones
    for (a.i in 1:m){
        ones <- which(zero.one[a.i, ] == 1)
        for (i in 1:m){
            if (any(zero.one[i, ones] == 1)){zero.one[i,ones] <- 1}
        } # end FOR i
    } # end FOR a.i

    ## build superclustering 1
    superclustering1 <- clustering1
    superclustering2 <- clustering2
    index <- rownames(zero.one)
    mapping1 <- list(); mapping2 <- list();
    i <- 1
    while (length(index) > 0){
        a.i <- index[1]
        cols <- which(zero.one[a.i, ] == 1)
        rows <- which(zero.one[, cols[1]] == 1)
        mapping2[[i]] <- names(cols); mapping1[[i]] <- names(rows)
        superclustering1[which(clustering1 %in% names(rows))] <-
            paste("S", i, sep = "")
        superclustering2[which(clustering2 %in% names(cols))] <-
            paste("T", i, sep = "")
        i <- i + 1
        remove <- which(index %in% names(rows))
        index <- index[-remove]
    } # end WHILE

    k <- length(mapping1)
    names(mapping1) <- paste("S", 1:k, sep = "")
    names(mapping2) <- paste("T", 1:k, sep = "")
    super.weights <- table(superclustering1, superclustering2)

    ## plotting
    if (plotting) {
        if (length(node.col) == 0) {node.col <- 1}
        if (length(edge.col) == 0) {edge.col <- 1}
        grav <- flatVSflat(superclustering1, superclustering2, 
            plotting = FALSE, h.min = h.min, max.iter = max.iter,
            greedy = FALSE)
        if (evenly) {
            coord1 <- (k:1)[order(grav$coord1, decreasing = TRUE)];
            coord2 <- (k:1)[order(grav$coord2, decreasing = TRUE)]
        } else {
            coord1 <- grav$coord1; coord2 <- grav$coord2
        }
        sc <- max(super.weights)
        super.weights.sc <- super.weights * line.wd / sc
        pmax <- max(c(coord1, coord2))
        pmin <- min(c(coord1, coord2))
        s <- (pmax - pmin) / 10

        if (horiz){
            coord1 <- (-coord1); coord2 <- (-coord2)
            pmax <- -pmax; pmin <- -pmin
            par(mar = c(1, 2, 1, 2) + 0.5)
            plot(x = c(coord1, coord2), y = c(rep(0, k), rep(1, k)),
                axes = FALSE, ty = "n",
                xlim = c(pmax - s, pmin + s),
                ylim = c(-3 * offset, 1 + 3 * offset), ...)
            myCex <- par("cex")

            # points
            points(x = coord1, y = rep(0, length(coord1)), cex = point.sz,
                pch = 19, col = node.col)
            points(x = coord2, y = rep(1, length(coord2)), cex = point.sz,
                pch = 19, col = node.col)

            # supercluster labels
            text(x = coord1, y = -offset,
                labels = rownames(super.weights), ...)
            text(x = coord2, y = 1 + offset,
                labels = colnames(super.weights), ...)
            par(cex = myCex * 0.7)

            for (sup in 1:k){
                # edges
                index <- which(super.weights.sc[sup, ] != 0)
                segments(coord1[sup], 0, coord2[index], 1,
                    lwd = super.weights.sc[sup, index], col = edge.col[1])

                # initial cluster labels
                text(x = coord1[sup], y = -2 * offset, labels =
                    paste( "(", paste(mapping1[[sup]], collapse = ","), ")",
                    sep = "" ), ...)
                text(x = coord2[sup], y = 1 + 2 * offset, labels =
                    paste( "(", paste(mapping2[[sup]], collapse = ","), ")",
                    sep="" ), ...)

                # sizes
                text(x = coord1[sup], y = -3 * offset, 
                    labels = paste( "size:", sum(super.weights[sup,]), 
                    sep = " " ), ...)
                text(x = coord2[sup], y = 1 + 3 * offset, 
                    labels = paste( "size:", sum(super.weights[,sup]),
                    sep=" " ), ...)
            } # end FOR sup

            par(cex = myCex)
        } else {  ## vertical
            par(mar = c(2, 1, 2, 1) + 0.5)
            plot(x = c(rep(0, k), rep(1, k)), y = c(coord1, coord2),
                axes = FALSE, xlim = c(-3 * offset, 1 + 3 * offset),
                ylim = c(pmin - s, pmax + s), ty = "n", ...)
            myCex <- par("cex")

            # points
            points(x = rep(0, length(coord1)), y = coord1, cex = point.sz,
                pch = 19, col = node.col)
            points(x = rep(1, length(coord2)), y = coord2, cex = point.sz,
                pch = 19, col = node.col)

            # supercluster labels
            text(x = -offset, y = coord1, labels = rownames(super.weights),
                ...)
            text(x = 1 + offset, y = coord2, labels = colnames(super.weights),
                ...)
            par(cex = myCex * 0.7)

            for (sup in 1:k){
                # edges
                index <- which(super.weights.sc[sup, ] != 0)
                segments(0, coord1[sup], 1, coord2[index],
                    lwd = super.weights.sc[sup, index], col = edge.col[1])

                # initial cluster labels
                text(x = -2 * offset, y = coord1[sup], labels = paste("(",
                    paste(mapping1[[sup]], collapse = ","), ")", sep = ""),
                    srt = 90, ...)
                text(x = 1 + 2 * offset, y = coord2[sup], labels = paste("(",
                    paste(mapping2[[sup]], collapse = ","), ")", sep = ""),
                    srt = 90, ...)

                # sizes
                text(x = -3 * offset, y = coord1[sup], 
                    labels = paste( "size:", sum(super.weights[sup,]), 
                    sep = " " ), ...)
                text(x = 1 + 3 * offset, y = coord2[sup], 
                    labels = paste( "size:", sum(super.weights[,sup]),
                    sep=" " ), ...)
            } # end FOR sup

            par(cex = myCex)
        } # end ELSE
    }  # end IF plotting

    return(list(s.clustering1 = superclustering1, s.clustering2 =
        superclustering2, merging1 = mapping1, merging2 = mapping2,
        weights = super.weights))
}
