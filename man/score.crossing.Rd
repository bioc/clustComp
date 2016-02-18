\name{score.crossing}
\alias{score.crossing}
\title{ Computation of the aesthetics-based score of the parent and the 
children trees }
\description{
\code{score.crossing} computes the value of the scoring function based on 
the aesthetics of the bi-graph formed when comparing a dendrogram and a 
flat clustering, for both the parent-tree and the children-tree; the 
children-tree consists of the same branches as the parent-tree, except for 
the parent node, that has been split and replaced by some of its 
descendants.}
\usage{
score.crossing(weight.1, weight.2, N.cross)
}
\arguments{
    \item{weight.1}{ a matrix of dimension \code{(m-1)xn} containing the
intersection sizes (edge weights) between branches in the parent-tree and 
clusters from the flat partitioning. The ordering of the rows and columns 
is irrelevant for the computation of the score.} 
    \item{weight.2}{ a matrix of dimension \code{(m+k)xn} containing the 
intersection sizes (edge weights) between branches in the children-tree and 
clusters from the flat partitioning. \code{k} takes on values in 0,1,...,L, 
where L is the maximum number of steps that the comparison algorithm is 
allowed to look ahead. The ordering of the rows corresponding to branches 
that are not descendants of the parent node must coincide with that of the 
matrix \code{weight.1} after discarding the parent node. The ordering of 
the columns is irrelevant for the computation of the score.} 
    \item{N.cross}{ the number of edge crossings induced in the subtree 
formed by the descendants of the parent-node.}
}
\value{
a list containing the following components:
    \item{sc1}{ the value of the scoring function for the parent-tree.}
    \item{sc2}{ the value of the scoring function for the children-tree.}
}
\details{ The decision to split a given parent-node is based on achieving a 
better score for the children-tree than for the parent-tree. In the case of 
\code{score.crossing}, a better score is reflected by a larger value of the 
scoring function, which rewards few thicker edges, penalises many smaller 
edges, and accounts for the number of edge crossings in the resulting 
bigraph. 
The descendants of the parent-node considered in the children-tree are its 
two children if no look-ahead is carried out; otherwise, the descendants 
will reach subsequent generations and their number will increase by one at 
each look-ahead step. }
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
\seealso{score.it, dyn.cross, flatVShier}
\examples{
    ### simulated data
    parent.clustering <- c(rep("B1", 5), rep("B2", 10), rep("B3", 10))
    # replace the branch 'B2' by its children 'B4' and 'B5'
    children.clustering <- c(rep("B1", 5), rep("B4", 3), rep("B5", 7),
        rep("B3", 10))
    flat.clustering <- c(rep(1, 6), rep(2, 6), rep(3, 4), rep(4, 9))
    # the ordering of flat clusters is '1','2','3' and '4'.
    parent.weights <- table(parent.clustering, flat.clustering)
    children.weights <- table(children.clustering, flat.clustering)
    descendant.cross <- dyn.cross(children.weights[c('B4', 'B5'), ])
    score.crossing(parent.weights, children.weights, descendant.cross) 
    ## better score for the parent.tree 
}
\keyword{ branch split }
