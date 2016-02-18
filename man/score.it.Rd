\name{score.it}
\alias{score.it}
\title{ Computation of the information theoretic-based score of the parent
and the children trees }
\description{
\code{score.it} computes the value of the scoring function based on 
information theory and the mutual information shared by a dendrogram and 
the flat clustering which is compared to, for both the parent-tree and the 
children-tree; the children-tree consists of the same branches as the 
parent-tree, except for the parent node, that has been split and replaced 
by some of its descendants.}
\usage{
score.it(weight.1, weight.2)
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
}
\value{
a list containing the following components:
    \item{sc1}{ the value of the scoring function for the parent-tree.}
    \item{sc2}{ the value of the scoring function for the children-tree.}
}
\details{ The decision to split a given parent-node is based on achieving a 
better score for the children-tree than for the parent-tree. In the case of 
\code{score.it}, a better score is reflected by a smaller value of the 
scoring function, which is related to the average length of the messages 
that encode the information about one clustering contained in the other. 
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
visualising relationships between hierarchical and flat gene expression 
data clusterings. \emph{Bioinformatics}, 21 (21), 3993-3999.
}
\seealso{score.crossing, flatVShier}
\examples{
    ### simulated data
    parent.clustering <- c(rep(1, 5), rep(2, 10), rep(3, 10))
    # replace the branch '2' by children '4' and '5'
    children.clustering<-c(rep(1,5),rep(4,3),rep(5,7),rep(3,10))
    flat.clustering <- c(rep(1, 6), rep(2, 6), rep(3, 4), rep(4, 9))
    score.it(table(parent.clustering, flat.clustering),
        table(children.clustering, flat.clustering)) 
    ## better score for the parent.tree 
}
\keyword{ branch split }
