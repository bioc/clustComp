\name{barycentre}
\alias{barycentre}
\title{ Computation of the barycentre-coordinate of a node connected to 
others in a bigraph}
\description{
\code{barycentre} provides an updated coordinate value for a node that is 
connected to nodes in a different layer of a bigraph. It is computed as the 
average of coordinates of the adjacent nodes, where weighted edges are 
considered as multi-edges collapsed into one.}
\usage{
barycentre(edge.weight, coordinates = NULL)
}
\arguments{
    \item{edge.weight}{ a vector containing the intersection sizes (edge 
weights) between a given node (from one of the partitionings) and all nodes 
in the other.} 
    \item{coordinates}{ a vector indicating the coordinates of the adjacent 
nodes. If it is not provided, then they are evenly spaced assuming a layout 
from top downwards.} 
}
\value{
    \item{position}{ a nsumber indicating the barycentre-coordinate of the 
node under consideration.  }
}
\details{ The node under consideration, from a given partitioning, is 
assigned a new position using the barycentre algorithm. The coordinates of 
the incident nodes are considered as many times as the corresponding edge 
weights indicate, and the new position for the node is given by the average 
over this set of coordinates. }
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
\seealso{flatVSflat, flatVShier}
\examples{
    # simulated data
    clustering1 <- c(rep(1, 5), rep(2, 10), rep(3, 10))
    clustering2 <- c(rep(1, 6), rep(2, 6), rep(3, 4), rep(4, 9))
    weights <- table(clustering1, clustering2)
    barycentre(weights[1, ], 1:4)
    barycentre(weights[1, ])
    barycentre(weights[2, ], 1:4)
    apply(weights, 1, barycentre, 1:4)
}
\keyword{ edge crossing }
