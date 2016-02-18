flatVShier <- function (tree, flat.clustering, flat.order = NULL,
    max.branches = 100, look.ahead = 2, pausing = TRUE, verbose = TRUE,
    h.min = 0.04, line.wd = 3, greedy = TRUE, greedy.colours = NULL,
    score.function = "crossing", expanded = FALSE, labels = NULL,
    cex.labels = 1, main = NULL) {
        
        
    ##### INITIAL SETUP  ######

    flat.clustering <- as.vector(flat.clustering)
    m <- 1    ## root of the tree
    n <- length(unique(flat.clustering))
    if (length(flat.order) == 0) {
        flat.order <- 1:n
        names(flat.order) <- unique(flat.clustering)
    }

    if (length(flat.order) != n) {stop("flat.order has an incorrect length...")}
    if ((is.character(score.function) == FALSE)|
        (!(score.function %in% c("crossing","it")))) {
            stop("score.function must be a character string,
            equal to 'crossing' or 'it'...")
    }
    
    no.genes <- length(flat.clustering)
    tree <- as.hclust(tree)
    tree.maxHeight <- max(tree$height)
    tree$height <- tree$height/tree.maxHeight
    if (max.branches > no.genes) {max.branches <- no.genes}

    ## initial tree: only one node (root)

    best.tree <- rep((no.genes - 1), no.genes) ## a vector stating the branch
        # to which each gene is assigned (the root at the first stage)
    best.order1 <- as.character(no.genes - 1) ## to use them as names/rownames
    best.order2 <- flat.order
    coord2.ifFinger <- list()
    best.coord1 <- 0.5  ## root node, only
    best.coord2 <- (1:n) - n / 2

    queue <- no.genes - 1    ## branches in queue to be split/explored
        #(by depth-first search)
    split <- c()     ##  to be a matrix storing by rows each explored
        # parent (col1) with its two children (cols 2 and 3)
    split.branch <- no.genes - 1   # root
    sub.branches <- tree$merge[split.branch, ] ## children of the split branch
        #(the two clusters agglomerated to form the split.branch)
    step <- no.genes - 1    # next branch to be analysed in the loop
    steps <- c()       # set of split branches

    subtree <- sub.branches # descendants of branch under consideration
        #to be split by looking ahead

    # initialize pointers for the look ahead part
    looking.ahead <- FALSE
    finger <- queue[1]
    finger.tree <- best.tree
    finger.distance <- 0; names(finger.distance) <- finger
    counter <- 0

    ###    start the iterations
    while ((length(queue) > 0) & (m < max.branches)) {

        ##########################
        ### update parent.tree ###
        ##########################

        parent.tree <- best.tree   ## tree
        parent.order1 <- best.order1  ## orders
        parent.order2 <- best.order2

        parent.coord1 <- best.coord1    ## coords
        parent.coord2 <- best.coord2

        weight.1 <- table(finger.tree, flat.clustering)
        split <- as.matrix(rbind(split, c(split.branch, sub.branches)))
        steps <- c(steps, step)

        ###########################
        ### build children.tree ###
        ###########################
        
        m <- m + 1   ## one more node in the children.tree after split
        children.tree <- parent.tree
        queue.1 <- sub.branches[1]
        queue.2 <- sub.branches[2]

        ## reassign genes in split branch to one of new children
        while (length(queue.1) > 0) {
            if (queue.1[1] < 0) {
                children.tree[abs(queue.1[1])] <- sub.branches[1]
                queue.1 <- queue.1[-1]
            } else {
                row <- queue.1[1]
                queue.1 <- insert(queue.1, position = 2, tree$merge[row,])
                queue.1 <- queue.1[-1]
            }  # end ELSE
        }  ## end WHILE queue1
        while ( length(queue.2) > 0) {
            if (queue.2[1] < 0) {
                children.tree[abs(queue.2[1])] <- sub.branches[2]
                queue.2 <- queue.2[-1]
            }  else {
                row <- queue.2[1]
                queue.2 <- insert(queue.2, 2, tree$merge[row,])
                queue.2 <- queue.2[-1]
            }  #  end ELSE
        } ## end WHILE queue2

        slot <- which(as.numeric(parent.order1) == split.branch)

        ## where to insert the new children in children.order1
        
        # create orders
        children.order1 <- parent.order1
        children.order1 <- insert(children.order1, slot, 0)
        children.order1 <- replace(children.order1, slot:(slot + 1),
            sub.branches)
        children.order2 <- parent.order2  ## (= best.order2)
            ## before barycentre
        children.coord2 <- parent.coord2
        m <- length(children.order1)
        children.coord1<- (1:m) - m / 2        # equally spaced
        weight.2 <- table(children.tree, flat.clustering)
            ## weights after split
        before.crossing <- dyn.cross(weight.2[children.order1,
            children.order2])
            ## before barycentre

        ###########################################
        #########  GRAVITATION ALGORITHM  #########
        ###########################################

        #################################
        ### barycentre; flat side (1) ###
        #################################

        ## update coords
        children.coord2 <- apply(weight.2[children.order1, ], 2,
            FUN = barycentre, coordinates = children.coord1)
            ### to update coordinates and order
        children.order2 <- order(children.coord2)
        # guarantee a distance at least h.min between nodes in flat layer
        for (j in 1:(n - 1)) {
            if (children.coord2[children.order2][j + 1] -
            children.coord2[children.order2][j] < h.min ) {
                children.coord2[children.order2][(j + 1):n] <-
                children.coord2[children.order2][(j + 1):n] + h.min
            } # end IF
        } # end FOR j

        current.crossing <- dyn.cross(weight.2[children.order1,children.order2])

        #######################
        ### swap; flat side ###
        #######################
        
        ## swap adjacent nodes in flat layer
        for (j in 1:(n - 1)){
            children.swapped.flat.order2 <- children.order2
            children.swapped.flat.order2[j:(j + 1)]<-
                children.swapped.flat.order2[(j + 1):j]
            swapped.flat.crossing <- dyn.cross(weight.2[children.order1,
                children.swapped.flat.order2])

            if (swapped.flat.crossing < current.crossing) {
                current.crossing <- swapped.flat.crossing
                children.order2 <- children.swapped.flat.order2
                children.coord2[children.swapped.flat.order2[j:(j + 1)]] <-
                    children.coord2[children.swapped.flat.order2[(j + 1):j]]
            } # end IF        
        } # end FOR

        if (current.crossing > before.crossing){
            best.crossing <- before.crossing
        } else {
            best.crossing <- current.crossing
            best.coord2 <- children.coord2
            best.order2 <- children.order2
        } # end ELSE

        #####################
        ### swap branches ###
        #####################

        children.swapped.tree.order1 <- children.order1
        children.swapped.tree.coord1 <- children.coord1
        index <- which( children.order1 %in% as.numeric(sub.branches) )
        children.swapped.tree.order1[index] <-
            children.swapped.tree.order1[rev(index)]
        before.crossing <- dyn.cross(weight.2[children.swapped.tree.order1,
            best.order2])

        #################################
        ### barycentre; flat side (2) ###
        #################################
        
        children.coord2 <- apply(weight.2[children.swapped.tree.order1,], 2,
            FUN=barycentre, coordinates = children.swapped.tree.coord1)
        children.order2 <- order(children.coord2)
        for (j in 1:(n - 1)) {
            if (children.coord2[children.order2][j + 1] -
                children.coord2[children.order2][j] < h.min ) {
                children.coord2[children.order2][(j + 1):n] <-
                    children.coord2[children.order2][(j + 1):n] + h.min
            } # end IF
        } # end FOR j
        current.crossing <- dyn.cross(weight.2[children.swapped.tree.order1,
            children.order2])

        #######################
        ### swap; flat side ###
        #######################

        for (j in 1:(n - 1)){
            children.swapped.flat.order2 <- children.order2
            children.swapped.flat.order2[j:(j + 1)] <-
                children.order2[(j + 1):j]
            swapped.flat.crossing <- dyn.cross(weight.2
                [children.swapped.tree.order1, children.swapped.flat.order2])
            if (swapped.flat.crossing < current.crossing) {
                current.crossing <- swapped.flat.crossing
                children.order2 <- children.swapped.flat.order2
                children.coord2[children.swapped.flat.order2[j:(j + 1)]] <-
                    children.coord2[children.swapped.flat.order2[(j + 1):j]]
            } # end IF        
        } # end FOR

        if (current.crossing > before.crossing) {
            best.swapped.crossing <- before.crossing
            swapped.order2 <- best.order2
            swapped.coord2 <- best.coord2
        } else  {
            best.swapped.crossing <- current.crossing
            swapped.order2 <- children.order2
            swapped.coord2 <- children.coord2
        } # end ELSE

        if (best.swapped.crossing < best.crossing){
            best.crossing <- best.swapped.crossing
            children.coord1 <- children.swapped.tree.coord1
            children.coord2 <- swapped.coord2
            children.order1 <- children.swapped.tree.order1
            children.order2 <- swapped.order2
            sub.branches[1:2] <- sub.branches[2:1]
            j <- which(subtree %in% sub.branches)
            subtree[j] <- subtree[rev(j)]
            split[nrow(split), 2:3] <- split[nrow(split), 3:2]
            
            ### reverse the original dendrogram for the expanded version
            tree$merge[split.branch,] <- rev(tree$merge[split.branch,])
            children.branches<-vector("list", 2)
            # compute clusters in tree
            for (cc in 1:2){
                v <- as.numeric(sub.branches[cc])
                if (v < 0) {
                    children.branches[[cc]]<- (-v)
                } else{
                    check <- tree$merge[v,]
    
                    while (length(check) > 0){
                        if (check[1] < 0) {
                            children.branches[[cc]] <-
                                c(children.branches[[cc]], -check[1])
                            check <- check[-1]
                        } else {
                            check <- insert(check, 2, tree$merge[check[1], ])
                            check <- check[-1]
                        } # end ELSE
                    } # end WHILE
                } # end ELSE
            } # end FOR
            
            group1<-sort(which( tree$order %in% children.branches[[1]] ))
            group2<-sort(which( tree$order %in% children.branches[[2]] ))
            lim1 <- min(c(group1, group2))
            lim2 <- max(c(group1, group2))
            if (lim1 > 1){
                if (lim2 < length(tree$order)){
                    tree$order <- tree$order[c(1:(lim1 - 1),group1,group2,
                        (lim2 + 1):length(tree$order))]
                } else {
                    tree$order <- tree$order[c(1:(lim1 - 1), group1, group2)]
                } # end ELSE
            } else {
                if (lim2 < length(tree$order)){
                    tree$order <- tree$order[c(group1, group2,
                        (lim2 + 1):length(tree$order))]
                } else {
                    tree$order <- tree$order[c(group1, group2)]
                } # end ELSE 
            } # end ELSE
        } else {
            children.coord2 <- best.coord2
            children.order2 <- best.order2
        } # end ELSE
            

        ################################################################
        ############          DRAW  THE  TREE       ####################
        ################################################################

        if (pausing) {
            dendrogram <- list(heights = tree$height[steps], branches = split)
            drawTreeGraph(weight.2, list(children.order1, children.order2),
                coordinates = list(children.coord1, children.coord2),
                dendrogram, line.wd = line.wd, main = paste("No.crossings = ",
                best.crossing, ", No.Branches = ", m), expanded = expanded,
                hclust.obj = tree, cex.labels = cex.labels, labels = labels)
        } # end IF pausing

        ####################################################################
        #########       DECISION ABOUT SPLITTING THE BRANCH      ###########
        ####################################################################

        ######## compute the score
        k <- which(children.order1 %in% subtree)
        subtree.order.1 <- children.order1[k]

        if (score.function == "crossing"){
            subtree.score <-score.crossing(weight.1, weight.2,
                dyn.cross(weight.2[subtree.order.1, children.order2]))
        } else {subtree.score <- lapply(score.it(weight.1, weight.2), "-")}
        
        if (looking.ahead) {score <- subtree.score$sc2 - finger.score.1}
        else {score <- subtree.score$sc2 - subtree.score$sc1}

        if (score >= 0) { # improvement

            ##############
            ### CASE A ###
            ##############
    
            # update queue
            if (queue[1] == (no.genes - 1)) {
                queue <- insert(queue, 2, tree$merge[step,])
            } else {queue <- insert(queue, 2, sub.branches)} # end ELSE

            if (verbose) {
                message(paste("Splitting the branch ", queue[1] ,
                    " in sub-branches ", queue[2] , " and ", queue[3],
                    " improves the score...", sep=""))
            } # end IF
            queue <- queue[-1]
            while (queue[1] < 0 & length(queue) > 0) {
                queue <- queue[-1]
            }  # remove leaves from the queue  
    
            if (looking.ahead){
                looking.ahead <- FALSE
                if (verbose) {"This split improves the scoring of
                    a previously not-allowed split..."}
                } # end IF looking-ahead
                if (length(queue) > 0) {
                    finger <- queue[1]
                    finger.tree <- children.tree
                    finger.distance <- 0; names(finger.distance) <- finger
                } # end IF update finger
                subtree <- tree$merge[finger, ]

                # update best
                best.tree <- children.tree
                best.order1 <- children.order1; best.order2 <- children.order2
                best.coord1 <- children.coord1; best.coord2 <- children.coord2
            }  else {  # no improvement
                if (pausing){
                    # mark the finger when no improvement
                    x <- tree$height[finger]
                    abline(v = -x, col = 2, lty = 2)
                } # end IF pausing
                # update counter
                counter <- counter + 1
                if (verbose) {message(paste("Splitting the branch ", queue[1],
                    " does not improve the score...", sep = "" ) )}
                if (looking.ahead == FALSE) {
                    # start a new look ahead process
                    looking.ahead <- TRUE
                    finger.score.1 <- subtree.score$sc1
                } # end IF new looking ahead
                # update generations
                finger.distance <- c(finger.distance, rep(counter, 2))
                names(finger.distance)[(length(finger.distance) - 1):
                    length(finger.distance)] <- tree$merge[queue[1], ]
                coord2.ifFinger[[length(coord2.ifFinger) + 1]] <- parent.coord2
                names(coord2.ifFinger) <- c(names(coord2.ifFinger)[1:
                    (length(coord2.ifFinger) - 1)], as.numeric(split.branch))

                if (counter <= look.ahead) {#further exploring unless leaves
                    if (verbose) {message("We look ahead one more step...")}

                    ##############
                    ### CASE B ###
                    ##############
                    
                    # case B.1 = two branches
                    if (all(sub.branches[1:2] > 0)){ # split branches normally
                        index <- which(subtree==sub.branches[1])
                        subtree <- insert(subtree, index + 1,
                            tree$merge[sub.branches[1], ])
                        subtree <- subtree[-index]
                    } else {
                        
                        if (all(sub.branches[1:2] < 0)){ # can't split branches
                        # even though the max number of steps ahead is not
                        # reached case B.2 = two leaves: undo the last split
                        
                        # update the subtree
                        index <- which(subtree %in% sub.branches)
                        subtree[index[1]] <- split.branch
                        subtree <- subtree[-(index[2])]
                        
                        # update the split/steps info
                        split <- as.matrix(rbind(c(), split[-nrow(split), ]))
                        steps <-steps[-length(steps)]
                        
                        # update children
                        children.tree <- parent.tree
                        children.coord1 <- parent.coord1
                        children.coord2 <- parent.coord2
                        children.order1 <- parent.order1
                        children.order2 <- parent.order2
                        m <- length(children.coord1)
                        A <- queue[1]
                        A.generation <- finger.distance[as.character(A)]
                        same.generation <- names(which(finger.distance ==
                            A.generation))
                        A.sibbling <- same.generation[which(same.generation %in%
                            queue[-1])]
                        index <- which(A.sibbling < 0)
                        if (length(index) > 0) {
                            A.sibbling <- A.sibbling[-index]
                        }
                        
                        if (length(A.sibbling) > 0) {
                            index <- which(subtree == A.sibbling[1])
                            subtree <- insert(subtree, index + 1,
                                tree$merge[as.numeric(A.sibbling[1]), ])
                            subtree <- subtree[-index]
                        } else {
                            
                            continue <- TRUE
                            
                            if (A == finger) {
                                continue <- FALSE
                                looking.ahead <- FALSE
                                coord2.ifFinger <- list()
                                
                                # update finger
                                index <- which(queue[-1] > 0) + 1
                                if (length(index) > 0) {
                                    finger <- queue[index[1]]
                                    finger.tree <- children.tree
                                    finger.distance <- 0
                                    names(finger.distance) <- finger
                                    subtree <- tree$merge[finger, ]
                                } # end IF
                                
                                if (verbose){message("No more splits allowed")}
                            } # end IF reach finger back
                            
                            while (continue){
                                A.parent <-which (tree$merge == A,
                                    arr.ind = TRUE)[, 1]
                                index <- which(subtree %in% tree$merge[A.parent,
                                    ])
                                subtree <- insert(subtree, index[1], A.parent)
                                subtree <- subtree[-(index + 1)]
                                A.p.generation <- finger.distance[as.character
                                    (A.parent)]
                                same.generation <- names(which
                                    (finger.distance == A.p.generation))
                                A.parent.sibbling <- same.generation[which
                                    (same.generation %in% queue[-1])]
                                index <- which(A.parent.sibbling < 0)
                                if (length(index) > 0) {A.parent.sibbling <-
                                    A.parent.sibbling[-index]}
                                if (length(A.parent.sibbling) > 0) {
                                    continue <- FALSE
                                    index <- which(subtree == A.parent.sibbling)
                                    subtree<-insert(subtree, index + 1,
                                        tree$merge
                                        [as.numeric(A.parent.sibbling), ])
                                    subtree <- subtree[-index]
                                } # end IF parent's sibbling 

                                A <- A.parent
                                #### back to the parent node
                                children.tree[children.tree %in% tree$merge
                                    [A.parent,]] <- A.parent
                                index <- which(children.order1 %in% tree$merge
                                    [A.parent, ])
                                children.order1[index[1]] <- A.parent
                                children.order1 <- children.order1[-(index[2])]
                                m <- m - 1
                                children.coord1 <- (1:m) - m / 2
                                index <- which(names(coord2.ifFinger) ==
                                    A.parent)
                                children.coord2 <- coord2.ifFinger[[index]]
                                children.order2 <- order(children.coord2)
                                index <- which(split[, 1] == A.parent)
                                split <- as.matrix(rbind(c(), split[-index,]))
                                steps <- steps[-index]
                
                                if (A == finger) {
                                    continue <- FALSE
                                    looking.ahead <- FALSE
                                    coord2.ifFinger <- list()
                                    # update finger
                                    index <- which(queue[-1] > 0) + 1
                                    if (length(index) > 0) {
                                        finger <- queue[index[1]]
                                        finger.tree <- children.tree
                                        finger.distance <- 0
                                        names(finger.distance) <- finger
                                        subtree <- tree$merge[finger, ]
                                    } # end IF                        
                                    
                                    if (verbose){message("No splits allowed.")}
                                } # end IF reach finger back
                            } # end WHILE continue
                        } # end ELSE    
                    } # end IF case B.2
                        
                    else { # case B.3 = one leaf, one branch
                        
                        if ((sub.branches[1] < 0) & (sub.branches[2] > 0)){
                            # can't split the leaf: will be removed later from
                            # the queue; split the branch...
                            index <- which(subtree == sub.branches[2])
                            subtree <- insert(subtree, index + 1,
                                tree$merge[sub.branches[2], ])
                            subtree <- subtree[-index]
                        } # end IF case B.3
                        
                        else{ # case B.4 = one branch, one leaf
                            index <- which(subtree == sub.branches[1])
                            subtree <- insert(subtree, index + 1,
                                tree$merge[sub.branches[1], ])
                            subtree <- subtree[-index]
                        } # end ELSE case B.4
                        
                    } # end ELSE                 
                } # end ELSE 
            
                if (queue[1] == (no.genes - 1)) {
                    queue <- insert(queue, 2, tree$merge[step, ])}
                else {queue <- insert(queue, 2, sub.branches)}  # end ELSE
                
                # update best
                best.tree <- children.tree
                best.order1 <- children.order1
                best.order2 <- children.order2
                best.coord1 <- children.coord1
                best.coord2 <- children.coord2
            }  ### end IF   case B
                
            else {  ## no further exploring is permitted...
                
                if (verbose) {
                    message("Can't explore this descendant any further.")
                } # end IF verbose

                ##############
                ### CASE C ###
                ##############
                
                ## similar to case B.2: we can't go further,
                #thus we step back and keep on exploring until possible...

                # update the subtree
                index <- which(subtree %in% sub.branches)
                subtree[index[1]] <- split.branch
                subtree <- subtree[-(index[2])]

                # update the split/steps info
                split <- as.matrix(rbind(c(), split[-nrow(split), ]))
                steps <- steps[-length(steps)]

                # update children
                children.tree <- parent.tree
                children.coord1 <- parent.coord1
                children.coord2 <- parent.coord2
                children.order1 <- parent.order1
                children.order2 <- parent.order2
                m <- length(children.coord1)
        
                A <- queue[1]

                A.generation <- finger.distance[as.character(A)]
                same.generation <- names(which(finger.distance == A.generation))
                A.sibbling <- same.generation[which(same.generation %in%
                    queue[-1])]
                index <- which(A.sibbling < 0)
                if (length(index) > 0){A.sibbling <- A.sibbling[-index]}
    
                if (length(A.sibbling) > 0) { # there's a explorable sibbling
                    index <- which(subtree == A.sibbling[1])
                    subtree <- insert(subtree, index + 1,
                        tree$merge[as.numeric(A.sibbling[1]), ])
                    subtree <- subtree[-index]
                } # end IF sibbling
                
                else { # there's no sibbling to explore
                    continue <- TRUE
                    if (A == finger) {
                        continue <- FALSE
                        looking.ahead <- FALSE
                        coord2.ifFinger <- list()
                        
                        # update finger
                        index <- which(queue[-1] > 0) + 1
                        if (length(index) > 0) {
                            finger <- queue[index[1]]
                            finger.tree <- children.tree
                            finger.distance <- 0
                            names(finger.distance) <- finger
                            subtree <- tree$merge[finger, ]
                        } # end IF
                        
                        if (verbose) {message("No more splits allowed...")}
                    } # end IF reach finger back
                    
                    while (continue){
                        A.parent <- which(tree$merge == A, arr.ind = TRUE)[, 1]
                        index <- which(subtree %in% tree$merge[A.parent, ])
                        subtree <- insert(subtree, index[1], A.parent)
                        subtree <- subtree[-(index + 1)]
                        A.p.generation <-
                            finger.distance[as.character(A.parent)]
                        same.generation <- names(which(finger.distance ==
                            A.p.generation))
                        A.parent.sibbling <- same.generation[which
                            (same.generation %in% queue[-1])]
                        index <- which(A.parent.sibbling < 0)
                        if (length(index) > 0){A.parent.sibbling <-
                            A.parent.sibbling[-index]}
                        if (length(A.parent.sibbling) > 0) {
                            continue <- FALSE
                            index <- which(subtree == A.parent.sibbling)
                            subtree <- insert(subtree, index + 1,
                                tree$merge[as.numeric(A.parent.sibbling), ])
                            subtree <- subtree[-index]
                        } # end IF parent's sibbling

                        A <- A.parent
                        
                        # back to the parent node
                        children.tree[children.tree %in%
                            tree$merge[A.parent, ]] <- A.parent
                        index <- which(children.order1 %in%
                            tree$merge[A.parent, ])
                        children.order1[index[1]] <- A.parent
                        children.order1 <- children.order1[-(index[2])]
                        m <- m - 1
                        children.coord1 <- (1:m) - m / 2
                        index <- which(names(coord2.ifFinger) == A.parent)
                        children.coord2 <- coord2.ifFinger[[index]]
                        children.order2 <- order(children.coord2)
                        index <- which(split[,1] == A.parent)
                        split <- as.matrix(rbind(c(), split[-index, ]))
                        steps <- steps[-index]
                        if (A == finger) {
                            continue <- FALSE
                            looking.ahead <- FALSE
                            coord2.ifFinger <- list()
                            if (verbose) {
                                message("No more splits are allowed.")}
                            
                            # update finger
                            index <- which(queue[-1] > 0) + 1
                            if (length(index) > 0) {
                                finger <- queue[index[1]]
                                finger.tree <- children.tree
                                finger.distance <- 0
                                names(finger.distance) <- finger
                                subtree <- tree$merge[finger, ]
                            } # end IF
                            
                        } # end IF reach finger back
                        
                    } # end WHILE continue
                    
                } # end ELSE no sibbling

                # update best
                best.tree <- children.tree
                best.order1 <- children.order1
                best.order2 <- children.order2
                best.coord1 <- children.coord1
                best.coord2 <- children.coord2
                
            } # end ELSE case C

            # update the queue for next loop
            queue <- queue[-1]
            while (queue[1] < 0 & length(queue) > 0) {queue <- queue[-1]}
                ## can't split leaves; remove them from the queue
            
        }  ### end ELSE   no improvement
            
        ### final updates for next loop
        step <- queue[1]
        split.branch <- queue[1]
        sub.branches <- tree$merge[step, ]
        counter <- finger.distance[as.character(queue[1])]
        names(counter) <- c()
        if (is.na(counter)){counter <- 0}
        if (pausing == TRUE) {pause()}
    } #### end WHILE 
    
    ################  final plot   #################################
    best.dendrogram <- list(heights = tree$height[steps], branches = split)
    weight.2 <- table(best.tree, flat.clustering)

    if (length(main) == 0) {
        main <- paste("Optimal cutoff - ", m, " branches", sep = "")
    }
    
    final.coords <- drawTreeGraph(weight.2, list(best.order1, best.order2),
        coordinates = list(best.coord1, best.coord2), best.dendrogram,
        main = main, dot = FALSE, line.wd = line.wd, expanded = expanded,
        hclust.obj = tree, cex.labels = cex.labels, labels = labels)

    ##################################################
    #################    GREEDY      #################
    ##################################################

    if (length(steps)==0) {greedy<-FALSE}

    if (greedy){
        Greedy <- SCmapping(best.tree, flat.clustering, plotting = FALSE)
        N <- length(unique(Greedy$s.clustering1)) # number of superclusters
        L <- length(greedy.colours) # number of colours provided
        if (L == 0) {greedy.colours <- 1:8; L <- 8}
        if (N > L * 5){
            print("Too many superclusters to be distinctly shown using the
                colours provided...")
        }
        else {
            
            if ((N > L)){
                greedy.colours <- rep(greedy.colours,ceiling(N/L))
            } #end IF
        
            greedy.symbols <- c(rep(21, L), rep(22, L), rep(23, L), rep(24, L),
                rep(25, L))[1 : N]
            for (p in 1:length(Greedy$merging1)){
                
                # label tree
                i <- which(best.order1 %in% Greedy$merging1[[p]])
                
                if (expanded) {
                    y <- final.coords$b.coord[i]
                    points(rep(final.coords$x.coords[1], length(y)), y,
                        col = greedy.colours[p], cex = 3,
                        pch = greedy.symbols[p])
                } # end IF
                
                else {
                    y <- best.coord1[i]
                    points(rep( final.coords$x.coords[1] ,length(y)), y,
                        col = greedy.colours[p], cex = 3,
                        pch = greedy.symbols[p])
                } # end ELSE
                                
                # label flat
                i <- which(colnames(weight.2) %in% Greedy$merging2[[p]])
                if (expanded){
                    y <- final.coords$f.coord[i]
                    points(rep(  final.coords$x.coords[2] ,length(y)), y,
                        col = greedy.colours[p], cex = 3,
                        pch = greedy.symbols[p])
                } # end IF
                
                else {
                    y <- best.coord2[i]
                    points(rep(final.coords$x.coords[2], length(y)), y,
                        col = greedy.colours[p], cex = 3,
                        pch = greedy.symbols[p])
                } # end ELSE
                
            } # end FOR
            
        } # end ELSE
        
        return(list(tree.partition = best.tree,
            tree.s.clustering = Greedy$s.clustering1,
            flat.s.clustering = Greedy$s.clustering2,
            tree.merging = Greedy$merging1,
            flat.merging = Greedy$merging2,
            dendrogram = tree))
    } # end IF
    
    else{
        return(list(tree.partition = best.tree, dendrogram = tree))
    } # end ELSE
} 