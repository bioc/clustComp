drawTreeGraph <- function(weight, current.order, coordinates, tree,
        dot = TRUE, line.wd = 3, main = NULL, expanded = FALSE,
        hclust.obj = NULL, flat.obj = NULL, labels = NULL, 
        cex.labels = 1, expression = NULL, layout = NULL, 
        ramp = NULL, bar1.col = NULL, bar2.col = NULL) {

    if (length(names(coordinates[[2]])) == 0) {
        names(coordinates[[2]]) <- colnames(weight)
    }

    heatmap <- (length(expression)>0)
    if (expanded == FALSE) {heatmap <- FALSE}
    else {if (length(hclust.obj)==0) 
        {stop('The full dendrogram must be provided to plot its expanded 
            version...')}}

    # determine layout (x axis) values for plot
    if (heatmap) {
        if (length(layout) < 3) {layout <- c(1,3,5)}
        else {layout <- sort(layout[1:3])}
    }
    else {
        if (length(layout) < 2) {layout <- c(1,3)}
        else {layout <- sort(layout[1:2])}
    }
    if (layout[1] <= 0) {stop("The layout values must be positive...")}

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

            B <- layout[length(layout)]
            space <- B - layout[length(layout)-1]
            A <- B - space/2

            plot( mySeq, ty = "n", axes = FALSE, ylab = '', xlab = '',
                main = main, xlim = c(-max(Height), B + space),
                ylim = c(min(mySeq), max(mySeq)) )
            myDots<-c()

            for (i in 1:(N-1)){
                if (all(Merge[i,] < 0)){
                    leaves <- which(Order %in% -(Merge[i,]))
                    myDots <- as.matrix(rbind(myDots, c(mySeq[leaves[1]],
                        mySeq[leaves[2]])))
                    # draw tree    
                    lines(x = c( 0, -Height[i], -Height[i], 0 ),
                        y = c( myDots[i,1], myDots[i,1], myDots[i,2],
                        myDots[i,2] ), col = 4 )
                } else {
                    if (all(Merge[i,] > 0)){
                        myDots <- as.matrix(rbind(myDots,
                            c(mean(myDots[Merge[i,1],]),
                            mean(myDots[Merge[i,2],]) ) ) )
                        # draw tree    
                        lines( x = c( -Height[Merge[i,1]], -Height[i],
                            -Height[i], -Height[Merge[i,2]] ),
                            y = c( myDots[i,1], myDots[i,1], myDots[i,2],
                            myDots[i,2] ), col = 4 )
                    } else {
                        if (Merge[i,1] < 0) {
                            leaves <- which(Order == -(Merge[i,1]))
                            myDots <- as.matrix(rbind(myDots, c(mySeq[leaves],
                                mean(myDots[Merge[i,2],]) ) ))
                            # draw tree    
                            lines(x = c( 0, -Height[i], -Height[i],
                                -Height[Merge[i,2]] ),
                            y = c( myDots[i,1], myDots[i,1], myDots[i,2],
                                myDots[i,2] ), col = 4 )
                        } else { # Merge[i,2]<0
                            leaves <- which(Order == -(Merge[i,2]))
                            myDots <- as.matrix(rbind(myDots,
                                c( mean(myDots[Merge[i,1],]), mySeq[leaves])))
                            # draw tree    
                            lines(x = c( -Height[Merge[i,1]], -Height[i],
                                -Height[i], 0),
                                y = c( myDots[i,1], myDots[i,1], myDots[i,2],
                                myDots[i,2] ), col = 4 )
                        } # end ELSE
                    } # end ELSE
                } # end ELSE 
            } # end FOR

            # write tree labels
            text(x = 0.1, y = mySeq, Labels[Order], adj = 0,
                cex = cex.labels)
            sep <- (mySeq[2] - mySeq[1]) / 2

            ## draw rectangles representing collapsed branches
            if (length(bar1.col)==0) {bar1.col <- 2}
            limits <- c(min(mySeq) - sep, max(mySeq) + sep )
            even.branch.coordinates <- (limits[2] + limits[1]) / 2
            rect(rep(layout[length(layout) - 1] + space/15, 2), limits[1], 
                rep(layout[length(layout) - 1] + space/5, 2), 
                limits[2], col = bar1.col[1], border = NA)

            ## draw heatmap            
            if (heatmap) {
                no.genes <- ncol(expression)
                no.samples <- nrow(expression)
                NN <- no.genes * no.samples
                if (length(ramp) == 2) { 
                    myColor <- colorRampPalette(c(ramp[1], ramp[2]))
                    } 
                else {myColor <- colorRampPalette(c('aliceblue', 
                    'darkcyan'))
                    }
                image(x = seq(layout[1], layout[2], length.out = no.genes + 1),
                    y = seq(min(limits), max(limits), 
                    length.out = no.samples + 1),
                    z = t(expression[Order,]), add = TRUE, col = myColor(NN))
                }

            # draw collapsed branches
            points(x = A, 
                y = even.branch.coordinates, pch = 19)

            # draw flat nodes
            points(x = rep(B, nodes.2), y = coordinates[[2]], pch = 19)

            # size of flat clusters
            if (length(bar2.col)==0) {bar2.col <- (1:nodes.2)}
            text(x = 0.2 + B, y = sort(coordinates[[2]]),
                labels = paste("size:", weight[current.order[[2]]]),
                cex = cex.labels * 0.9, adj = 0)

           # labels
            text(x = rep(B + space/2, nodes.2), y = coordinates[[2]], 
                labels = paste("F(", colnames(weight)[sort(current.order[[2]])],
                ")", sep = ''),
                #labels = colnames(weight)[sort(current.order[[2]])], 
                adj = 0, cex = cex.labels, col = bar2.col)

            ## draw rectangles representing flat clusters
            bar2.colours <- match(flat.obj, sort(unique(colnames(weight))))
            rect(rep( B + space, 
                nodes.2 + 1), 
                c(min(mySeq-sep),sep + mySeq[1:N-1]), 
                rep(B + 5*space/3, 
                nodes.2 + 1), 
                sep+mySeq, 
                col = (bar2.col[bar2.colours])[Order],
                border = NA)

            # draw edges
            sc <- max(weight); weight.sc <- weight * line.wd / sc
            segments(rep(A, nodes.2), rep(even.branch.coordinates, nodes.2),
                rep(B, nodes.2), coordinates[[2]], lwd = weight.sc)
            if (dot) {
                bb <- tree$branches[nrow(tree$branches), 1]
                points(x = -Height[bb], y = mean(myDots[bb,]), col = 3,
                    cex = 1.1)
            } # end IF

            even.flat.coordinates <- coordinates[[2]]

        } else{  #############    ONE BRANCH, COLLAPSED
            B <- layout[length(layout)] 
            space <- B - layout[length(layout)-1]
            A <- B - space/2            

            plot(rep(B, nodes.2), coordinates[[2]], xlab = "", ylab = "",
                xlim = c(A - space, B + space), 
                ylim = c(min(min(coordinates[[1]],
                min(coordinates[[2]]))) - 0.1, max(max(coordinates[[1]],
                max(coordinates[[2]]))) + 0.1),
                main = main, xaxt = "n", yaxt = "n", pch = 19, axes = FALSE)

            # draw root
            points(x = A, y = 0.5, pch = 19)
            sc <- max(weight); weight.sc <- weight * line.wd / sc
            segments( rep(A, length(coordinates[[2]])), 0.5, 
                B, coordinates[[2]], lwd = weight.sc)

            # label flat clusters
            text(B + space/2,
                sort(coordinates[[2]]), labels = paste("F(",
                names(coordinates[[2]])[order(coordinates[[2]])],
                ")", sep = ""), cex = cex.labels, col = 2, adj = 1)

            # size of tree
            text(A - 0.2, 0.5, labels = paste("size:", 
                sum(weight)), cex = cex.labels * 0.9, adj = 1)

            # size of flat clusters
            text(B + 0.2, sort(coordinates[[2]]),
                labels = paste("size:", weight[current.order[[2]]]),
                cex = cex.labels * 0.9, adj = 0)
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
            B <- layout[length(layout)]
            space <- B - layout[length(layout)-1]
            A <- B - space/2

            plot(mySeq, ty = "n", axes = FALSE, ylab = '', xlab = '',
                xlim = c(-max(Height), B + space), 
                ylim = c(min(mySeq), max(mySeq)), main = main)
            myDots <- c()

            for (i in 1:(N-1)){
                if (all(Merge[i,] < 0)){
                    leaves <- which(Order %in% -(Merge[i,]))
                    myDots <- as.matrix(rbind(myDots, c(mySeq[leaves[1]],
                        mySeq[leaves[2]])))
                    # draw tree    
                    lines(x = c(0, -Height[i], -Height[i], 0),
                        y = c(myDots[i,1], myDots[i,1], myDots[i,2],
                        myDots[i,2] ), col = 4 )
                } else {
                    if (all(Merge[i,] > 0)){
                        myDots <- as.matrix(rbind(myDots,
                            c(mean(myDots[Merge[i,1],]),
                            mean(myDots[Merge[i,2],]) ) ) )
                        # draw tree    
                        lines(x = c(-Height[Merge[i,1]], -Height[i], -Height[i],
                            -Height[Merge[i,2]] ), y = c( myDots[i,1],
                            myDots[i,1], myDots[i,2], myDots[i,2] ), col = 4)
                    } else {
                        if (Merge[i,1] < 0) {
                            leaves <- which(Order == -(Merge[i,1]))
                            myDots <- as.matrix(rbind(myDots, c(mySeq[leaves],
                                mean(myDots[Merge[i,2],]) ) ) )
                            # draw tree
                            lines(x = c(0, -Height[i], -Height[i],
                                -Height[Merge[i,2]]), y = c(myDots[i,1],
                                myDots[i,1], myDots[i,2], myDots[i,2]), col = 4)
                        } else {  # Merge[i,2]<0
                            leaves <- which(Order == -(Merge[i,2]))
                            myDots <- as.matrix(rbind(myDots,
                                c(mean(myDots[Merge[i,1],]), mySeq[leaves])))
                            # draw tree    
                            lines(x = c( -Height[Merge[i,1]], -Height[i],
                                -Height[i],0), y = c(myDots[i,1], myDots[i,1],
                                myDots[i,2], myDots[i,2] ), col = 4 )
                        } # end ELSE
                    } # end ELSE
                } # end ELSE 
            } # end FOR

            # write tree labels
            text(x = 0.1, y = mySeq, Labels[Order], adj = 0,
                cex = cex.labels)
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

            ## draw heatmap

            if (heatmap) {
                no.genes <- ncol(expression)
                no.samples <- nrow(expression)
                NN <- no.genes * no.samples
                if (length(ramp) == 2) { 
                    myColor <- colorRampPalette(c(ramp[1], ramp[2]))
                    } 
                else {myColor <- colorRampPalette(c('aliceblue', 
                    'darkcyan'))
                    }
                image(x = seq(layout[1], layout[2], length.out = no.genes + 1),
                    y = seq(min(limits), max(limits), 
                    length.out = no.samples + 1),
                    z = t(expression[Order,]), add = TRUE, col = myColor(NN))
                }

            ## draw rectangles representing collapsed branches
            if (length(bar1.col)==0) {bar1.col <- 1:nodes.1}
            limits <- sort(limits)
            even.branch.coordinates <- (limits[2:(nodes.1 + 1)] +
                limits[1:nodes.1]) / 2
            rect(rep(layout[length(layout) - 1] + space/15, nodes.1 + 1), 
                limits[1:nodes.1], rep(layout[length(layout)-1] + space/5, 
                nodes.1 + 1), limits[2:(nodes.1 + 1)], col = bar1.col, 
                border = NA)

            # draw collapsed branches
            points(x = rep(A, nodes.1),
                y = even.branch.coordinates, pch = 19)

            # draw flat nodes
            even.flat.coordinates <- (min(mySeq) + h.flat *
                seq(1, nodes.2 * 2, by = 2))[order(current.order[[2]])]
            points(x = rep(B, nodes.2),
                y = even.flat.coordinates, pch = 19)

            # write size box2
            if (length(bar2.col)==0) {bar2.col <- (1:nodes.2)}
            text(x = 0.2 + B, 
                y = sort(even.flat.coordinates),
                labels = paste("size:", apply(weight[, current.order[[2]]],
                2, sum)), cex = cex.labels * 0.9, adj = 0)

            # labels
            text(x = rep(B + space/2, 
                nodes.2), y = even.flat.coordinates, 
                labels = paste("F(", colnames(weight)[sort(current.order[[2]])],
                ")", sep = ''), 
                #labels = colnames(weight)[sort(current.order[[2]])],
                adj = 0,
                cex = cex.labels, col = bar2.col)

            ## draw rectangles representing flat clusters
            bar2.colours <- match(flat.obj, sort(unique(colnames(weight))))
            rect(rep( B + space, 
                nodes.2 + 1), 
                c(min(mySeq-sep),sep + mySeq[1:N-1]), 
                rep(B + 5*space/3, 
                nodes.2 + 1), 
                sep+mySeq, 
                col = (bar2.col[bar2.colours])[Order],
                border = NA)

            # draw edges
            sc <- max(weight)
            weight.sc <- weight * line.wd / sc

           for (i in 1:nodes.1){
                index <- which(weight.sc[current.order[[1]][i], ] != 0)
                segments(A, even.branch.coordinates[i],
                    B, even.flat.coordinates[index],
                    lwd = weight.sc[current.order[[1]][i], index])
            } # end FOR
            if (dot) {
                bb <- tree$branches[nrow(tree$branches), 1]
                points(x = -Height[bb], y = mean(myDots[bb,]),
                col = 3, cex = 1.1)
            } # end IF

        } else{    ############ collapsed tree
            B <- layout[length(layout)] 
            space <- B - layout[length(layout)-1]
            A <- B - space/2            

            # plot the branches
            plot(rep(A, nodes.1), coordinates[[1]],
                xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                xlim = c(-1, B + space),
                ylim = c(min(min(coordinates[[1]], min(coordinates[[2]]))) -
                0.1, max(max(coordinates[[1]], max(coordinates[[2]]))) + 0.1 ),
                main = main, pch = 19, axes = FALSE)

            # plot the flat clusters
            points(rep(B, nodes.2), coordinates[[2]], pch = 19)

           # branch labels
            text(0.1, sort(coordinates[[1]]), labels = paste( "B(",
                current.order[[1]], ")", sep = ""), cex = cex.labels,
                col = 2, adj = 0)

            # label flat clusters
            text(B + space/2,
                sort(coordinates[[2]]), labels = paste( "F(",
                names(coordinates[[2]])[order(coordinates[[2]])], ")",
                sep = ""), cex = cex.labels, col = 2, adj = 0)

            # write size box1
            text(A - 0.2, sort(coordinates[[1]]),
            labels = paste("size:", apply(weight[current.order[[1]], ], 1,
                sum)), cex = cex.labels * 0.9, adj = 1)

            # write size box2
            text(B + 0.2, sort(coordinates[[2]]),
                labels = paste("size:", apply(weight[, current.order[[2]]],
                2, sum)), cex = cex.labels * 0.9, adj = 0)

            # draw edges
            sc <- max(weight); weight.sc <- weight * line.wd / sc
            for (i in 1:nodes.1){
                index <- which(weight.sc[current.order[[1]][i], ] != 0)
                segments(A, sort(coordinates[[1]])[i], 
                    B, coordinates[[2]][index],
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
            x.coords = c(A, B)))
    } else {
        return(list(b.coord = sort(coordinates[[1]]),
            f.coord = coordinates[[2]],
            x.coords = c(A, B)))
    } # end ELSE
}

