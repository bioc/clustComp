flatVSflat <- function(flat1, flat2, coord1 = NULL, coord2 = NULL, 
    max.iter = 24, h.min = 0.1, plotting = TRUE, horiz = FALSE, 
    offset = 0.1, line.wd = 3, point.sz = 2, evenly = FALSE, 
    greedy = TRUE, greedy.colours = NULL, main = "", xlab = "", 
    ylab = "", col = NULL, ...){

    weights <- as.matrix(table(flat1, flat2))
    m <- nrow(weights)
    n <- ncol(weights)
    if (length(rownames(weights)) == 0) {
        rownames(weights) <- paste("A", 1:m, sep = "")
    }
    if (length(colnames(weights)) == 0) {
        colnames(weights) <- paste("B", 1:n, sep = "")
    }

    if (length(coord1) == 0) {
        coord1 <- m:1 - m / 2
        names(coord1) <- rownames(weights)
    } else {coord1 <- as.vector(coord1)}
    if (length(coord2) == 0) {
        coord2 <- n:1 - n / 2
        names(coord2) <- colnames(weights)
    } else {coord2 <- as.vector(coord2)}

    if ((any(names(coord1) != rownames(weights))) | (any(names(coord2) !=
        colnames(weights)))) {
            stop("The coordinates and the contingency table labels do not
                match...")
    }

    current.order1 <- order(coord1, decreasing = TRUE)
    current.order2 <- order(coord2, decreasing = TRUE)
    current.crossings <- initial.crossings <-
        dyn.cross(weights[current.order1, current.order2])

    current.coord1 <- coord1
    current.coord2 <- coord2
    continue <- TRUE
    it.counter <- 0
    if (it.counter >= max.iter) {continue <- FALSE}
    if ((m==1)|(n==1)){
        continue <- FALSE
        current.crossings <- 0
    } # end IF

    while (continue){
        previous.order1 <- current.order1
        previous.order2 <- current.order2

        ### find the barycentre for each node in each layer
        ### layer 1
        new.coord1 <- apply(weights, 1, FUN = barycentre,
            coordinates = current.coord2)
        new.order1 <- order(new.coord1, decreasing = TRUE)

        ### correction of distances if too close
        for (i in 1:(m-1)) {
            if (new.coord1[new.order1[i]] - new.coord1[new.order1[i + 1]] <
            h.min) {
                new.coord1[new.order1[(i+1):m]] <-
                    new.coord1[new.order1[(i+1):m]] - h.min
            } # end IF
        } # end FOR

        ### compute the number of crossings; update if there is any improvement
        new.crossings <- dyn.cross(weights[new.order1, current.order2])
        improved <- (new.crossings < current.crossings)
        if (improved) {
            current.order1 <- new.order1
            current.coord1 <- new.coord1
            current.crossings <- new.crossings
        } # end IF improved

        # swap consecutive nodes
        improved <- TRUE
        while (improved){
            improved <- FALSE
            for (i in 1:(m - 1)){
                new.order1 <- current.order1
                new.order1[i:(i + 1)] <- new.order1[(i + 1):i]
                new.crossings <- dyn.cross(weights[new.order1, current.order2])
                if (new.crossings < current.crossings) {
                    current.crossings <- new.crossings
                    current.order1 <- new.order1
                    swapping <- c(names(current.coord1[ current.order1[i] ]),
                        names(current.coord1[ current.order1[i + 1] ]) )
                    current.coord1[swapping] <- current.coord1[rev(swapping)]
                    improved <- TRUE
                } # end IF
            } # end FOR
        } # end WHILE

        ### layer 2
        new.coord2 <- apply(weights, 2, FUN = barycentre,
            coordinates = current.coord1)
        new.order2 <- order(new.coord2, decreasing = TRUE)

        ### correction of distances if too close
        for (j in 1:(n - 1)) {
            if (new.coord2[new.order2[j]] -
                new.coord2[new.order2[j + 1]] < h.min){
                new.coord2[new.order2[(j + 1):n]] <-
                    new.coord2[new.order2[(j + 1):n]] - h.min
            } # end IF
        } # end FOR j

        ### compute the number of crossings; update if there is any improvement
        new.crossings <- dyn.cross(weights[current.order1, new.order2])
        improved <- (new.crossings < current.crossings)
        if (improved) {
            current.order2 <- new.order2
            current.coord2 <- new.coord2
            current.crossings <- new.crossings
        } # end IF

        # swap consecutive nodes
        improved <- TRUE
        while (improved){
            improved <- FALSE
            for (j in 1:(n - 1)){
                new.order2 <- current.order2
                new.order2[j:(j + 1)] <- new.order2[(j + 1):j]
                new.crossings <- dyn.cross(weights[current.order1, new.order2])
                if (new.crossings < current.crossings) {
                    current.crossings <- new.crossings
                    current.order2 <- new.order2
                    swapping <- c( names(current.coord2[ current.order2[j] ]), 
                        names( current.coord2[ current.order2[j + 1] ]) )
                    current.coord2[swapping] <- current.coord2[rev(swapping)]
                    improved <- TRUE
                } # end IF
            } # end FOR
        } # end WHILE

        it.counter <- it.counter + 1
        if (all(previous.order1 == current.order1) &
            all(previous.order2 == current.order2)) {continue <- FALSE}
        if (it.counter >= max.iter) {continue <- FALSE}
    } # end WHILE

    if (plotting){
        if (evenly) {
            current.coord1 <- c(m:1)[order(current.order1)] - m / 2
            current.coord2 <- c(n:1)[order(current.order2)] - n / 2
        } # end IF evenly

        pmax <- max(c(current.coord1, current.coord2))
        pmin <- min(c(current.coord1, current.coord2))
        s <- (pmax - pmin) / 10
        par(mar = c(1, 1, 1, 1) + 0.5)
        if (length(col) == 0) {col <- 1}
        if (horiz){
            current.coord1 <- (-current.coord1)
            current.coord2 <- (-current.coord2)
            pmax <- -pmax
            pmin <- -pmin
            plot(x = c(current.coord1, current.coord2),
                y=c(rep(0, length(current.coord1)),
                rep(1, length(current.coord2))),
                axes = FALSE, main = main, xlab = xlab, ylab = ylab,
                xlim = c(pmax - s, pmin + s),
                ylim = c(-2 * offset, 1 + 2 * offset), ty="n")

            # edges
            sc <- max(weights); weights.sc <- weights*line.wd / sc
            for (i in 1:m){
                order <- which(weights.sc[i, ] != 0)
                segments(current.coord1[i], 0, current.coord2[order], 1,
                    lwd = weights.sc[i, order], col = col)
            } # end FOR

            # points
            points(x = current.coord1, y = rep(0, length(current.coord1)),
                cex = point.sz, pch = 19, col = col)
            points(x = current.coord2, y = rep(1, length(current.coord2)),
                cex = point.sz, pch = 19, col = col)

            # labels
            text(x = current.coord1, y = -offset, labels = rownames(weights),
                col = col, ...)
            text(x = current.coord2, y = 1 + offset, labels = colnames(weights),
                col = col, ...)

            # greedy
            if (greedy){
                Greedy <- SCmapping(flat1, flat2, plotting = FALSE)
                # number of superclusters
                N <- length(unique(Greedy$s.clustering1)) 
                # number of colours provided
                L <- length(greedy.colours) 
                if (L == 0) {greedy.colours <- 1:8; L <- 8}
                if (N > L * 5){
                    print("Too many superclusters to be distinctly shown using
                        the colours provided...")
                } # end IF
                else {
                    if ((N > L)){
                        greedy.colours <- rep(greedy.colours,ceiling(N/L))
                    } #end IF

                    greedy.symbols <- c(rep(21, L), rep(22, L), 
                        rep(23, L), rep(24, L), rep(25, L))[1 : N]
                    for (p in 1:length(Greedy$merging1)){

                        # label flat1
                        i <- which(rownames(weights) %in% Greedy$merging1[[p]])                        
                        points(x = current.coord1[i], 
                            y = rep(0, length(current.coord1[i])),
                            col = greedy.colours[p], cex = point.sz + 2, 
                            pch = greedy.symbols[p])

                        # label flat2
                        i <- which(colnames(weights) %in% Greedy$merging2[[p]])
                        y <- current.coord2[i]
                        points(x = current.coord2[i], 
                            y = rep(1, length(current.coord2[i])),
                            col = greedy.colours[p], cex = point.sz + 2, 
                            pch = greedy.symbols[p])
                    }   #  end FOR p
                } # end ELSE
            } # end IF greedy 
        } else {
            plot(x = c(rep(0, m), rep(1, n)),
                y = c(current.coord1, current.coord2),
                axes = FALSE, main = main, xlab = xlab, ylab = ylab,
                xlim = c(-2 * offset, 1 + 2 * offset),
                ylim = c(pmin - s, pmax + s), ty = "n")

            # edges
            sc <- max(weights); weights.sc <- weights * line.wd / sc
            for (i in 1:m){
                index <- which(weights.sc[i, ] != 0)
                segments(0, current.coord1[i], 1, current.coord2[index],
                    lwd = weights.sc[i, index], col = col)
            } # end FOR

            # points
            points(x = rep(0, m), y = current.coord1, cex = point.sz, pch = 19,
                col = col)
            points(x = rep(1, n), y = current.coord2, cex = point.sz, pch = 19,
                col = col)

            # labels
            text(x = -offset, y = current.coord1, labels = rownames(weights),
                col = col, ...)
            text(x = 1 + offset, y = current.coord2, labels = colnames(weights),
                col = col, ...)

            # greedy
            if (greedy){
                Greedy <- SCmapping(flat1, flat2, plotting = FALSE)
                # number of superclusters
                N <- length(unique(Greedy$s.clustering1)) 
                # number of colours provided
                L <- length(greedy.colours) 
                if (L == 0) {greedy.colours <- 1:8; L <- 8}
                if (N > L * 5){
                    print("Too many superclusters to be distinctly shown using
                        the colours provided...")
                } # end IF
                else {
                    if ((N > L)){
                        greedy.colours <- rep(greedy.colours,ceiling(N/L))
                    } #end IF

                    greedy.symbols <- c(rep(21, L), rep(22, L), 
                        rep(23, L), rep(24, L), rep(25, L))[1 : N]
                    for (p in 1:length(Greedy$merging1)){

                        # label flat1
                        i <- which(rownames(weights) %in% Greedy$merging1[[p]])                        
                        points(x = rep(0, length(current.coord1[i])),
                            y = current.coord1[i], 
                            col = greedy.colours[p], cex = point.sz + 2, 
                            pch = greedy.symbols[p])

                        # label flat2
                        i <- which(colnames(weights) %in% Greedy$merging2[[p]])
                        y <- current.coord2[i]
                        points(x = rep(1, length(current.coord2[i])),
                            y = current.coord2[i],
                            col = greedy.colours[p], cex = point.sz + 2, 
                            pch = greedy.symbols[p])
                    }   #  end FOR p

                } # end ELSE
            } # end IF greedy 

        } # end ELSE

    } # end IF plotting

    if (greedy) {
        return(list(icross = initial.crossings, fcross = current.crossings,
            coord1 = current.coord1, coord2 = current.coord2, 
            s.clustering1 = Greedy$s.clustering1,
            s.clustering2 = Greedy$s.clustering2,
            merging1 = Greedy$merging1,
            merging2 = Greedy$merging2))
    } else {
        return(list(icross = initial.crossings, fcross = current.crossings,
            coord1 = current.coord1, coord2 = current.coord2))
    }
}
