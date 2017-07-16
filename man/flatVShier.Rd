\name{flatVShier}
\alias{flatVShier}
\title{ Comparison of a hierarchical and a flat clusterings}
\description{
    \code{flatVShier} carries out the comparison and visualisation of the 
relationships between a hierarchical and a flat clusterings. The 
hierarchical one is shown either as a complete or pruned tree, whose 
collapsed branches are nodes on the left hand side layer of a bi-graph. The 
flat clusters are represented on the right hand side. Branches and flat 
clusters are connected with edges, whose thickness represents the number of 
elements common to both sets. The number of edge crossings is minimised 
using the barycentre algorithm on the right hand side; also, the children 
corresponding to the last split in the dendrogram when exploring it by 
depth-first search are swapped if this decreases the number of crossings.}
\usage{
flatVShier(tree, flat.clustering, flat.order = NULL, look.ahead = 2, 
    score.function = "crossing", expanded = FALSE, expression = NULL, 
    greedy = TRUE, layout = NULL, pausing = TRUE, verbose = TRUE, 
    greedy.colours = NULL, h.min = 0.04, labels = NULL, max.branches = 100,
    line.wd = 3, cex.labels = 1, main = NULL, ramp = NULL, bar1.col = NULL, 
    bar2.col = NULL)
}
\arguments{
    \item{tree}{ an \code{hclust} object, or structure that can be converted 
to hclust object, corresponding to a data set of size N.}
    \item{flat.clustering}{ a vector of length N containing the labels for 
each of the N objects clustered.}
    \item{flat.order}{ an optional vector containing the initial ordering 
for the flat clusters (from bottom upwards).}
    \item{look.ahead}{ the number of steps allowed to look further after 
finding a parent-node whose score is better that that of its children.}
    \item{score.function}{ a string specifying whether the decision to 
split a given branch is based on the aesthetics ('crossing') or on the 
information theory ('it').}
    \item{expanded}{ a Boolean parameter to state whether the hierarchical 
tree is displayed complete or pruned. FALSE by default.}
    \item{expression}{ a matrix containing the expression data from which
the dendrogram and the flat partioning were obtained. If provided, and 
argument \code{expanded} is set to TRUE, the heatmap of the data is 
represented.} 
    \item{greedy}{ a Boolean argument; if TRUE, the branches produced by 
the optimal cutoffs are used to construct superclusters that will be mapped 
to superclusters on the flat side with the greedy algorithm in 
\code{SCmapping}.} 
    \item{layout}{ a vector containing 3 or 2 components, depending on
whether the heatmap is plotted or not. Each coordinate provides the value
on the X axis for locating the heatmap, the colour bar representing pruned
branches, and the bigraph. Provided vectors are increasingly ordered to 
keep the layout of the plot. By default, the X values are 1, 3, (5).}
    \item{pausing}{ a Boolean argument; if TRUE, each step in the 
comparison is plotted and followed by a pause.}
    \item{verbose}{ a Boolean argument; if TRUE, the situation in each 
iteration is described.}
    \item{greedy.colours}{ an optional vector containing the colours for 
each of the superclusters if the greedy algorithm is used; if the length of 
this vector is smaller than that of the number of resulting superclusters 
p, then it is recycled.}
    \item{h.min}{ minimum separation between nodes in the flat layer; if 
the barycentre algorithm sets two nodes to be less than this distance 
apart, then the second node and the following ones are shifted (upwards).}
    \item{labels}{ an optional vector containing labels for the leaves of 
the tree.}
    \item{max.branches}{ an integer stating the maximum number of branches 
allowed in the dendogram.}
    \item{line.wd}{  a numerical parameter that fixes the width of the 
thickest edge(s); the rest are drawn proportionally to their weights; 3 by 
default.} 
    \item{cex.labels}{ a number indicating the magnification for the labels 
of the flat clusters and the branches in the hierarchical tree.}
    \item{main}{ an optional character string for the title of the plot.}
    \item{ramp}{ a vector with two components containing the two colours
used to define the palette for the heatmap.}
    \item{bar1.col}{ a vector of integers or character strings indicating 
the colours to be used in the hierarchical coloured bar, drawn on the left
hand side, if the dendrogram is fully expanded.}
    \item{bar2.col}{ a vector of integers or character strings indicating 
the colours to be used in the flat coloured bar, drawn on the right hand
side, if the dendrogram is fully expanded.}
}
\value{
    a list of components including:
    \item{tree.partition}{ a vector of length N stating the branch each 
element belongs to.}
    \item{tree.s.clustering}{ a vector of length N stating the supercluster 
on the tree side each element belongs to. If greedy=FALSE this component is 
not returned.}
    \item{flat.s.clustering}{ a vector of length N stating the supercluster 
on the flat side each element belongs to. If greedy=FALSE this component is 
not returned.}
    \item{tree.merging}{ a list of p components; the j-th element contains 
the labels of the tree that have been merged to produce the j-th 
supercluster. If greedy=FALSE this component is not returned.}
    \item{flat.merging}{ a list of p components; the j-th element contains 
the labels of the flat clusters that have been merged to produce the j-th 
supercluster. If greedy=FALSE this component is not returned.}
    \item{dendrogram}{ an hclust object with the appropriate ordering of 
the branches to minimise the number of crossings.}
}
\details{ 
The method cuts different branches of the tree at 'optimal' levels, which 
may be different at different branches, to find the best matches with the 
flat clustering. The method explores the tree depth-first, starting from 
the root. In each iteration the goal is to decide whether the branch under 
consideration (the parent node) is to be split. To that end, the user 
selects a scoring function based on the bi-graph aesthetics ('crossing') or 
an alternative based on mutual information between the flat and 
hierarchical clusterings ('it'). The selected score is first computed for 
the parent-node. Next it is replaced by its children; the barycentre 
algorithm, with the swapping strategy, is used on the flat side of the bi-
graph. Later, the children in the tree are swapped and the positions on the 
flat side are likewise updated. The best score obtained by any of these 
layouts in the children tree is compared to the score of the parent-tree, 
sp. If it is better, then the splitting is allowed and the tree is 
subsequently explored. Otherwise, the splitting is discarded, unless it is 
allowed to look ahead. In that case, the score for the tree with one of the 
children of the parent-node replaced by its own children is compared to sp; 
this is repeated until we get a better score for the children or until the 
maximum number of looking-ahead steps is reached.
After the optimal cut-offs are found, it is possible to run a greedy 
algorithm to determine sets of clusters from each side which have a large 
overlap. These sets, referred to as superclusters, determine the mapping 
between the two clusterings.
}
\author{
Aurora Torrente \email{aurora@ebi.ac.uk}
and
Alvis Brazma \email{brazma@ebi.ac.uk}
}
\references{
Torrente, A. \emph{et al.} (2005). A new algorithm for comparing and 
visualizing relationships between hierarchical and flat gene expression 
data clusterings. \emph{Bioinformatics}, 21 (21), 3993-3999.
}
\seealso{flatVSflat, barycentre, score.crossing, score.it, SCmapping}
\examples{
    # simulated data
    set.seed(0)
    dataset <- rbind(matrix(rnorm(20), 5, 4), sweep(matrix(rnorm(24), 6, 4),
        2, 1:4, "+"))
    tree <- hclust(dist(dataset))
    # two clusters
    flat <- kmeans(dataset,2)$cluster
    collapsed1 <- flatVShier(tree, flat, pausing = FALSE)
    # four clusters
    flat<-kmeans(dataset, 4)$cluster
    collapsed2 <- flatVShier(tree, flat)

    ## expanded tree
    # no heatmap
    expanded1 <- flatVShier(tree, flat, pausing = FALSE, score.function = "it",
        expanded = TRUE, bar1.col = c("red", "blue", "cyan", "magenta"),
        bar2.col = c("gray", "yellow", "orange", "purple"))
    # with heatmap
    expanded2 <- flatVShier(tree, flat, pausing = FALSE, score.function = "it",
        expanded = TRUE, expression = dataset,
        bar1.col = c("red", "blue", "cyan", "magenta"),
        bar2.col = c("gray", "yellow", "orange", "purple"))    
}
\keyword{ clustering comparison }
