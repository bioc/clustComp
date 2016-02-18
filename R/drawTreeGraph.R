drawTreeGraph <- function(weight, current.order, coordinates, tree,
        dot = TRUE, line.wd = 3, main = NULL, expanded = FALSE,
        hclust.obj = NULL, labels = NULL, cex.labels = 1) {
            
    if (length(names(coordinates[[2]])) == 0) {
        names(coordinates[[2]]) <- colnames(weight)
    }
    
    if (length(tree$heights) == 0) {
        message("No branch should be split...")
        nodes.2 <- length(current.order[[2]])
        
        if (expanded){   ############    ONE BRANCH, EXPANDED
            Merge <- hclust.obj$merge
            Height <- hclust.obj$height
            max.height <- max(Height)
            Height <- Height/max.height
            Order <- hclust.obj$order
            N <- length(Order)
        
            if (length(labels) > 0){
                Labels <- labels
                if (length(Labels) != length(Order)) {
                    stop("Incorrect number of labels...")
                }
            } else {
                Labels <- hclust.obj$labels
                if (length(Labels) == 0) {
                    Labels <- 1:N
                }
            }
            
            mySeq <- seq(min(coordinates[[2]]), max(coordinates[[2]]),
                length.out = N)
            space <- max(nchar(Labels)) * 0.05 * cex.labels
            space.flat <- max(nchar(current.order[[2]])) * 0.05 * cex.labels
            space.size1 <- nchar(as.character(N)) * 0.05 * cex.labels
            space.size2 <- max(nchar(sum(weight[current.order[[2]]]))) *
                0.05 * cex.labels
            plot( mySeq, ty = "n", axes = FALSE, ylab = '', xlab = '',
                main = main, xlim = c(-max(Height), 2 + space + space.flat +
                space.size1 + space.size2) )
            myDots<-c()
            for (i in 1:(N-1)){
                if (all(Merge[i,] < 0)){
                    leaves <- which(Order %in% -(Merge[i,]))
                    myDots <- as.matrix(rbind(myDots, c(mySeq[leaves[1]],
                        mySeq[leaves[2]])))
                    lines(x = c( 0, -Height[i], -Height[i], 0 ),
                        y = c( myDots[i,1], myDots[i,1], myDots[i,2],
                        myDots[i,2] ), col = 4 )
                } else {
                    if (all(Merge[i,] > 0)){
                        myDots <- as.matrix(rbind(myDots,
                            c(mean(myDots[Merge[i,1],]),
                            mean(myDots[Merge[i,2],]) ) ) )
                        lines( x = c( -Height[Merge[i,1]], -Height[i],
                            -Height[i], -Height[Merge[i,2]] ),
                            y = c( myDots[i,1], myDots[i,1], myDots[i,2],
                            myDots[i,2] ), col = 4 )
                    } else {
                        if (Merge[i,1] < 0) {
                            leaves <- which(Order == -(Merge[i,1]))
                            myDots <- as.matrix(rbind(myDots, c(mySeq[leaves],
                                mean(myDots[Merge[i,2],]) ) ))
                            lines(x = c( 0, -Height[i], -Height[i],
                                -Height[Merge[i,2]] ),
                            y = c( myDots[i,1], myDots[i,1], myDots[i,2],
                                myDots[i,2] ), col = 4 )
                        } else {  # Merge[i,2]<0
                            leaves <- which(Order == -(Merge[i,2]))
                            myDots <- as.matrix(rbind(myDots,
                                c( mean(myDots[Merge[i,1],]), mySeq[leaves])))
                            lines(x = c( -Height[Merge[i,1]], -Height[i],
                                -Height[i], 0),
                                y=c( myDots[i,1], myDots[i,1], myDots[i,2],
                                myDots[i,2] ), col = 4 )
                        } # end ELSE
                    } # end ELSE
                } # end ELSE 
            } # end FOR
            text(x = 0.1, y = mySeq, Labels[Order], adj = 0,
                cex = cex.labels * 0.7)
            sep <- (mySeq[2] - mySeq[1]) / 2

            ## draw rectangles representing collapsed branches
            limits <- c(min(mySeq) - sep, max(mySeq) + sep )
            even.branch.coordinates <- (limits[2] + limits[1]) / 2
            rect(rep(0.15 + space, 2), limits[1], rep(0.25 + space, 2),
                limits[2], col = 2, border = NA)

            # draw collapsed branches
            points(x = 0.35 + space + space.size1, y = even.branch.coordinates,
                pch = 19)

            # draw flat nodes
            points(x = rep(1.3 + space + space.size1, nodes.2),
                y = coordinates[[2]], pch = 19)

            # size of flat clusters
            text(1.4 + space + space.size1, sort(coordinates[[2]]),
                labels = paste("size:", weight[current.order[[2]]]),
                cex = cex.labels * 0.7, adj = 0, col = 2)

            # labels
            text(x = rep(1.7 + space + space.size1 + space.flat + space.size2,
                nodes.2), y = coordinates[[2]], labels = paste("F(",
                sort(current.order[[2]]), ")", sep = '' ), adj = 0,
                cex = cex.labels)

            # draw edges
            sc <- max(weight); weight.sc <- weight * line.wd / sc
            segments(rep(0.35 + space + space.size1, nodes.2),
                rep(even.branch.coordinates, nodes.2),
                rep(1.3 + space + space.size1, nodes.2),
                coordinates[[2]], lwd = weight.sc)
            if (dot) {
                bb <- tree$branches[nrow(tree$branches), 1]
                points(x = -Height[bb], y = mean(myDots[bb,]), col = 3,
                    cex = 1.1)
            } # end IF
            
            even.flat.coordinates <- coordinates[[2]]

        } else{       #############    ONE BRANCH, COLLAPSED
            space <- max(nchar(current.order[[1]])) * 0.05 * cex.labels
            space.flat <- max(nchar(current.order[[2]])) * 0.05 * cex.labels
            space.size1 <- nchar(as.character(length(hclust.obj$order))) *
                0.05 * cex.labels
            space.size2 <- max(nchar(weight)) * 0.05 * cex.labels
            plot(rep( 1.3 + space.size1, nodes.2), coordinates[[2]],
                xlab = "", ylab = "",
                xlim = c(0.25, 2 + space.flat + space.size1 + space.size2),
                ylim = c(min(min(coordinates[[1]],
                min(coordinates[[2]]))) - 0.1, max(max(coordinates[[1]],
                max(coordinates[[2]]))) + 0.1),
                main = main, xaxt = "n", yaxt = "n", pch = 19, axes = FALSE)

            # draw root
            points(x = 0.35 + space.size1, y = 0.5, pch = 19)
            sc <- max(weight); weight.sc <- weight * line.wd / sc
            segments( rep(0.35 + space.size1, length(coordinates[[2]])), 0.5,
                1.3 + space.size1, coordinates[[2]], lwd = weight.sc)
                
            # label flat clusters
            text(1.7 + space.size1 + space.flat + space.size2,
                sort(coordinates[[2]]), labels = paste("F(",
                names(coordinates[[2]])[order(coordinates[[2]])],
                ")", sep = ""), cex = cex.labels, col = 2, adj = 1)

            # size of tree
            text(0.25 + space.size1,  0.5, labels = paste("size:", sum(weight)),
                cex = cex.labels * 0.7, adj = 1)

            # size of flat clusters
            text(1.4 + space.size1, sort(coordinates[[2]]),
                labels = paste("size:", weight[current.order[[2]]]),
                cex = cex.labels * 0.7, adj = 0)
        } # end ELSE
        
    } else{   ###################     MORE THAN ONE BRANCH
        tree.maxHeight <- max(tree$height)
        tree.heights <- tree$heights/tree.maxHeight
        tree.branch <- tree$branch
        current.order[[1]] <- as.character(current.order[[1]])
        nodes.1 <- length(current.order[[1]])
        nodes.2 <- length(current.order[[2]])
        minimum <- min(c(coordinates[[1]], coordinates[[2]]))
        maximum <- max(c(coordinates[[1]], coordinates[[2]]))
        new.split <- tree.branch[nrow(tree.branch), 1]

        if (expanded) {   #########    EXPANDED
            Merge <- hclust.obj$merge
            Height <- hclust.obj$height
            max.height <- max(Height)
            Height <- Height / max.height
            Order <- hclust.obj$order
            N <- length(Order)

            if (length(labels) > 0){
                Labels <- labels
                if (length(Labels) != N) {stop("Incorrect number of labels.")}
            } else {
                Labels <- hclust.obj$labels
                if (length(Labels) == 0) {Labels <- 1:N}
            }
            
            mySeq <- seq(1, nodes.1, length.out = N) - nodes.1 / 2
            h.flat <- (max(mySeq) - min(mySeq)) / (nodes.2 * 2)
            space <- max(nchar(Labels)) * 0.05 * cex.labels
            space.flat <- max(nchar(current.order[[2]])) * 0.05 * cex.labels
            space.size1 <- max(nchar(apply(weight[current.order[[1]],], 1,
                sum))) * 0.05 * cex.labels
            space.size2 <- max(nchar(apply(weight[, current.order[[2]]],
                2,sum))) * 0.05 * cex.labels
            plot(mySeq, ty = "n", axes = FALSE, ylab = '', xlab = '',
                xlim = c(-max(Height), 2 + space + space.flat + space.size1 +
                space.size2), main = main)
            myDots <- c()

            for (i in 1:(N-1)){
                if (all(Merge[i,] < 0)){
                    leaves <- which(Order %in% -(Merge[i,]))
                    myDots <- as.matrix(rbind(myDots, c(mySeq[leaves[1]],
                        mySeq[leaves[2]])))
                    lines(x = c(0, -Height[i], -Height[i], 0),
                        y = c(myDots[i,1], myDots[i,1], myDots[i,2],
                        myDots[i,2] ), col = 4 )
                } else {
                    if (all(Merge[i,] > 0)){
                        myDots <- as.matrix(rbind(myDots,
                            c(mean(myDots[Merge[i,1],]),
                            mean(myDots[Merge[i,2],]) ) ) )
                        lines(x = c(-Height[Merge[i,1]], -Height[i], -Height[i],
                            -Height[Merge[i,2]] ), y = c( myDots[i,1],
                            myDots[i,1], myDots[i,2], myDots[i,2] ), col = 4)
                    } else {
                        if (Merge[i,1] < 0) {
                            leaves <- which(Order == -(Merge[i,1]))
                            myDots <- as.matrix(rbind(myDots, c(mySeq[leaves],
                                mean(myDots[Merge[i,2],]) ) ) )
                            lines(x = c(0, -Height[i], -Height[i],
                                -Height[Merge[i,2]]), y = c(myDots[i,1],
                                myDots[i,1], myDots[i,2], myDots[i,2]), col = 4)
                        } else {  # Merge[i,2]<0
                            leaves <- which(Order == -(Merge[i,2]))
                            myDots <- as.matrix(rbind(myDots,
                                c(mean(myDots[Merge[i,1],]), mySeq[leaves])))
                            lines(x = c( -Height[Merge[i,1]], -Height[i],
                                -Height[i],0), y = c(myDots[i,1], myDots[i,1],
                                myDots[i,2], myDots[i,2] ), col = 4 )
                        } # end ELSE
                    } # end ELSE
                } # end ELSE 
            } # end FOR
            
            text(x = 0.1, y = mySeq, Labels[Order], adj = 0,
                cex = cex.labels * 0.7)
            sep <- (mySeq[2] - mySeq[1]) / 2
            limits <- min(mySeq) - sep
            clusters <- vector("list", nodes.1)
            
            # compute clusters in tree and draw red dots
            for (cc in 1:nodes.1){
                v <- as.numeric(current.order[[1]][cc])
                if (v < 0) {
                    clusters[[cc]] <- (-v)
                    points(x = 0, y = mySeq[max(which(Order == (-v)))],
                        pch = 19, col = 2)
                } else {
                    points( x = -Height[v], y = mean(myDots[v,]), pch = 19,
                        col = 2)
                    check <- Merge[v,]
                    while (length(check) > 0){
                        if (check[1] < 0) {
                            clusters[[cc]] <- c(clusters[[cc]], -check[1])
                            check <- check[-1]
                        } else {
                            check <- insert(check, 2, Merge[check[1],])
                            check <- check[-1]
                        } # end ELSE
                    } # end WHILE
                } # end ELSE
                
                # compute limits for rectangles
                limits <- c( limits, sep +
                    mySeq[max(which(Order %in% clusters[[cc]]))] )
            } # end FOR
            
            ## draw rectangles representing collapsed branches
            limits <- sort(limits)
            even.branch.coordinates <- (limits[2:(nodes.1 + 1)] +
                limits[1:nodes.1]) / 2
            rect(rep(0.15 + space, nodes.1 + 1), limits[1:nodes.1],
                rep(0.25 + space, nodes.1 + 1), limits[2:(nodes.1 + 1)],
                col = 1:nodes.1, border = NA)
                
            # draw collapsed branches
            points(x = rep(0.35 + space + space.size1, nodes.1),
                y = even.branch.coordinates, pch = 19)
                
            # draw flat nodes
            even.flat.coordinates <- (min(mySeq) + h.flat *
                seq(1, nodes.2 * 2, by = 2))[order(current.order[[2]])]
            points(x = rep(1.3 + space + space.size1, nodes.2),
                y = even.flat.coordinates, pch = 19)
                
            # write size box2
            text(x = 1.4 + space + space.size1, y = sort(even.flat.coordinates),
                labels = paste("size:", apply(weight[, current.order[[2]]],
                2, sum)), cex = cex.labels * 0.7, adj = 0)
                
            # labels
            text(x = rep(1.7 + space + space.size1 + space.size2, nodes.2),
                y = even.flat.coordinates, labels = paste("F(",
                sort(current.order[[2]]), ")", sep = ''), adj = 0,
                cex = cex.labels, col = 2)
                
            # draw edges
            sc <- max(weight)
            weight.sc <- weight * line.wd / sc
            for (i in 1:nodes.1){
                index <- which(weight.sc[current.order[[1]][i], ] != 0)
                segments(0.35 + space + space.size1, even.branch.coordinates[i],
                    1.3 + space + space.size1, even.flat.coordinates[index],
                    lwd = weight.sc[current.order[[1]][i], index])
            } # end FOR
            
            if (dot) {
                bb <- tree$branches[nrow(tree$branches), 1]
                points(x = -Height[bb], y = mean(myDots[bb,]),
                col = 3, cex = 1.1)
            } # end IF
            
        } else{    ############ collapsed tree
            space <- max(nchar(current.order[[1]])) * 0.05 * cex.labels
            space.flat <- max(nchar(current.order[[2]])) * 0.05 * cex.labels
            space.size1 <- max(nchar(apply(weight[current.order[[1]], ], 1,
                sum))) * 0.05 * cex.labels + 0.2
            space.size2 <- max(nchar(apply(weight[, current.order[[2]]],
                2, sum))) * 0.05 * cex.labels
                
            # plot the branches
            plot(rep(0.35 + space + space.size1, nodes.1), coordinates[[1]],
                xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                xlim = c(-1, 2 + space + space.flat + space.size1 +
                space.size2),
                ylim = c(min(min(coordinates[[1]], min(coordinates[[2]]))) -
                0.1, max(max(coordinates[[1]], max(coordinates[[2]]))) + 0.1 ),
                main = main, pch = 19, axes = FALSE)

            # plot the flat clusters
            points(rep(1.3 + space + space.size1, nodes.2), coordinates[[2]],
                pch = 19)
            
            # branch labels
            text(0.1, sort(coordinates[[1]]), labels = paste( "B(",
                current.order[[1]], ")", sep = ""), cex = cex.labels,
                col = 2, adj = 0)
            # label flat clusters
            text(1.7 + space + space.size1 + space.flat + space.size2,
                sort(coordinates[[2]]), labels = paste( "F(",
                names(coordinates[[2]])[order(coordinates[[2]])], ")",
                sep = ""), cex = cex.labels, col = 2, adj = 0)
                
            # write size box1
            text(0.25 + space + space.size1, sort(coordinates[[1]]),
            labels = paste("size:", apply(weight[current.order[[1]], ], 1,
                sum)), cex = cex.labels * 0.7, adj = 1)
                
            # write size box2
            text(1.4 + space + space.size1, sort(coordinates[[2]]),
                labels = paste("size:", apply(weight[, current.order[[2]]],
                2, sum)), cex = cex.labels * 0.7, adj = 0)
                
            # draw edges
            sc <- max(weight); weight.sc <- weight * line.wd / sc
            for (i in 1:nodes.1){
                index <- which(weight.sc[current.order[[1]][i], ] != 0)
                segments(0.35 + space + space.size1, sort(coordinates[[1]])[i],
                    1.3 + space + space.size1, coordinates[[2]][index],
                    lwd = weight.sc[current.order[[1]][i], index])
            } # end FOR
            
            no.steps <- length(tree.heights)
            tree.coordinates <- cbind(current.order[[1]], rep(0, nodes.1),
                coordinates[[1]])
            tree.coordinates <- matrix(as.numeric(tree.coordinates),
                nrow(tree.coordinates), 3)
            for (step in no.steps:1) {
                y1 <- tree.coordinates[tree.coordinates[, 1] ==
                    tree.branch[step, 2], 3]
                y2 <- tree.coordinates[tree.coordinates[, 1] ==
                    tree.branch[step, 3], 3]

                # draw the dendrogram
                lines(x = c(-tree.coordinates[tree.coordinates[, 1] ==
                    tree.branch[step, 2], 2], -tree.heights[step],
                    -tree.heights[step],
                    -tree.coordinates[tree.coordinates[, 1] ==
                    tree.branch[step, 3], 2]),
                    y = c(y1, y1, y2, y2), col = 4)
                tree.coordinates[tree.coordinates[,1] ==
                    tree.branch[step,2],] <-
                c(tree.branch[step, 1], tree.heights[step], mean(c(y1, y2)))
                index <- which(tree.coordinates[, 1] == tree.branch[step,3],
                    arr.ind = TRUE)[1]
                tree.coordinates <- rbind(c(), tree.coordinates[-index, ])
                
                # label the branch that has been split with a green dot
                if (dot) {
                    index <- which(tree.coordinates[, 1] == new.split)
                    if (length(which(tree.coordinates[,1] == new.split)) != 0) {
                        points(x = -tree.coordinates[index, 2],
                        y = tree.coordinates[index, 3], col = 3, cex = 1.2)
                    } # end IF
                } # end IF dot
            } # end FOR step
        } # end ELSE
    } # end ELSE

    if (expanded) {
        return(list(b.coord = even.branch.coordinates,
            f.coord = even.flat.coordinates,
            x.coords = c(0.35 + space + space.size1,
            1.3 + space + space.size1)))
    } else {
        return(list(b.coord = sort(coordinates[[1]]),
            f.coord = coordinates[[2]],
            x.coords = c(0.35 + space + space.size1,
            1.3 + space + space.size1)))
    } # end ELSE
}