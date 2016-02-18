\name{dyn.cross}
\alias{dyn.cross}
\title{ Computation of the number of crossings in the bi-graph in a 
particular layout }
\description{
\code{dyn.cross} dynamically computes the number of edge crossings for a 
given node ordering of both layers of the bi-graph. This ordering is given 
by the rownames and colnames of the matrix of weights provided as argument. 
The use of this function is mainly internal.}
\usage{
dyn.cross(weights)
}
\arguments{
        \item{weights}{ a matrix containing the intersection sizes (edge 
weights) between clusters (nodes) from the first partitioning (left or top 
layer of the bi-graph) and the second partitioning (right or bottom layer 
of the bi-graph). The rownames and colnames of the matrix provide the 
ordering of the nodes in the first and second layer, respectively.} 
}
\value{
        \item{crossings}{ the number of weighted crossings in the layout 
provided as argument.  }
}
\details{The number of crossings in the weighted bigraph is computed by 
considering multiedges between connected nodes. The specific implementation 
uses a modification, adapted to the case of having multiple edges, of the 
dynamic programming algorithm developed by Nagamochi and Yamada for 
counting the number of crossings in non-weighted graphs, which has reduced 
computational cost. }
\author{
Aurora Torrente \email{aurora@ebi.ac.uk}
and
Alvis Brazma \email{brazma@ebi.ac.uk}
}
\references{
Nagamochi H. and Yamada N. (2004) Counting edge crossings in a 2-layered 
drawing. \emph{Information Processing Letters}. 91, 221-225.
}
\seealso{barycentre, flatVSflat, flatVShier}
\examples{
# simulated data
    clustering1 <- c(rep(1, 5), rep(2, 10), rep(3, 10))
    clustering2 <- c(rep(1, 6), rep(2, 6), rep(3, 4), rep(4, 9))
    dyn.cross(table(clustering1, clustering2)) # no crossings
    dyn.cross(table(clustering1, clustering2)[c(2, 1, 3), ]) 
}
\keyword{ edge crossing }
