\name{flatVSflat}
\alias{flatVSflat}
\title{ Comparison of two flat clusterings}
\description{
\code{flatVSflat} carries out the comparison and visualisation of two flat 
clusterings. The nodes in each partitioning are represented as nodes in the 
two layers of a bi-graph. The sizes of the intersection between clusters 
are reflected in the edge thickness. The number of edge crossings is 
minimised heuristically using the barycentre algorithm alternatively on 
each side.}
\usage{
flatVSflat(flat1, flat2, coord1 = NULL, coord2 = NULL, max.iter = 24, 
    h.min = 0.1, plotting = TRUE, horiz = FALSE, offset = 0.1, line.wd = 3, 
    point.sz = 2, evenly = FALSE, greedy = TRUE, greedy.colours = NULL, 
    main = "", xlab = "", ylab = "", col = NULL, ...)
}
\arguments{
    \item{flat1}{ a vector of length N containing the labels for each of the 
N genes clustered according to the first flat method.} 
    \item{flat2}{ a vector of length N containing the labels for each of the 
N genes clustered according to the second flat method.} 
    \item{coord1}{ a vector indicating the coordinates of the nodes in the 
first layer of the bi-graph. If not provided, then the nodes are initially 
equally spaced.}
    \item{coord2}{ a vector indicating the coordinates of the nodes in the 
second layer of the bi-graph. If not provided, then the nodes are initially 
equally spaced.} 
    \item{max.iter}{an integer stating the maximum number of runs of the 
barycentre heuristic on both layers of the bi-graph.}
    \item{h.min}{ minimum separation between nodes in the same layer; if the 
barycentre algorithm sets two nodes to be less than this distance apart, 
then the second node and the following ones are shifted (downwards, in the 
vertical layout, and to the right, in the horizontal layout).}
    \item{plotting}{ a Boolean parameter which yields the bi-graph if TRUE.}
    \item{horiz}{ a Boolean argument for displaying a vertical (default) or 
horizontal layout.}
    \item{offset}{ a numerical parameter that sets the separation between the 
nodes and their labels. It is set to 0.1 by default.}
    \item{line.wd}{ a numerical parameter that fixes the width of the thickest 
edge(s); the rest are drawn proportionally to their weights; 3 by default.}
    \item{point.sz}{ a numerical parameter that fixes the size of the nodes in 
the bigraph; 2 by default.}
    \item{evenly}{ a Boolean parameter; if TRUE the coordinate values are 
ignored, and the nodes are drawn evenly spaced, according to the ordering 
obtained by the algorithm. It is set to FALSE by default.}
    \item{greedy}{ a Boolean argument; if set to TRUE the greedy algorithm is
used to construct a one-to-one mapping between superclusters on both sides.}
    \item{greedy.colours}{ a vector of integers or character strings defining 
the colours that will be used to visualise the superclusters. If not provided 
colours are selected from the default palette of size 8. If the length of 
this vector is smaller than that of the number of resulting superclusters 
p, then it is recycled, and different symbols are used.}
    \item{main}{ graphical parameter as in \code{plot}.}
    \item{xlab}{ graphical parameter as in \code{plot}.}
    \item{ylab}{ graphical parameter as in \code{plot}.}
    \item{col}{ graphical parameter as in \code{plot}.}
    \item{...}{ further graphical parameters.}
}
\value{
a list of components including:
    \item{icross}{ the number of edge crossings before running the barycentre 
algorithm.}
    \item{fcross}{ the number of edge crossings after running the barycentre 
algorithm.}
    \item{coord1}{ a vector containing the coordinates for each node in the 
first layer.}
    \item{coord2}{ a vector containing the coordinates for each node in the 
second layer.}
    \item{s.clustering1}{ a vector of length N stating the supercluster 
in the first layer each element belongs to. If greedy=FALSE this component is 
not returned.}
    \item{s.clustering2}{ a vector of length N stating the supercluster 
in the second layer each element belongs to. If greedy=FALSE this component 
is not returned.}
    \item{merging1}{ a list of p components; the j-th element contains the 
labels from the first clustering that have been merged to produce the j-th 
supercluster. If greedy=FALSE this component is not returned.}
    \item{merging2}{ a list of p components; the j-th element contains the 
labels from the second clustering that have been merged to produce the j-th 
supercluster. If greedy=FALSE this component is not returned.}
}
\details{ 
As the iterations of the algorithm run the coordinates of the nodes in a 
single layer are updated. For a given partition, each node is assigned a 
new position, the gravity-centre, using the barycentre algorithm; then, the 
nodes in the corresponding layer are reordered according to the new 
positions. If the gravity-centres cause two consecutive nodes to be less 
than \code{h.min} apart, the coordinates of the second and all the 
following ones are shifted.
Additionally, to improve the results of the algorithm the following 
strategy is also used after running the barycentre algorithm on each side: 
consecutive nodes are swapped if this transposition leads to a reduction in 
the number of edge crossings.
The algorithm runs until there is no improvement in the number of crossings 
or until the maximum number of iterations is reached.
The rownames and colnames of matrix \code{weights} contain the cluster 
labels.
The ordering in the layout is over-imposed by the coordinate values, 
therefore, the names (in the coordinates) and row-/col-names (in the 
contingency table) should coincide.
}
\author{
Aurora Torrente \email{aurora@ebi.ac.uk}
and
Alvis Brazma \email{brazma@ebi.ac.uk}
}
\references{
Eades, P.  \emph{et al.} (1986). On an edge crossing problem. Proc. of 9th 
Australian Computer Science Conference, pp. 327-334.

Gansner, E.R.  \emph{et al.} (1993). A technique for drawing directed 
graphs. IEEE Trans. on Software Engineering, 19 (3), 214-230.

Garey, M.R.  \emph{et al.} (1983). Crossing number in NP complete. SIAM J. 
Algebraic Discrete Methods, 4, 312-316.

Torrente, A. \emph{et al.} (2005). A new algorithm for comparing and 
visualizing relationships between hierarchical and flat gene expression 
data clusterings. \emph{Bioinformatics}, 21 (21), 3993-3999.
}
\seealso{flatVShier, barycentre}
\examples{
    ### simulated data; two superclusters
    clustering1 <- c(rep(1, 5), rep(2, 10), rep(3, 10))
    clustering2 <- c(rep(1, 4), rep(2,15), rep(1, 6))
    # two colours
    flatVSflat(clustering1, clustering2, greedy.colours = c("red","blue"))
    # one colour, two symbols
    flatVSflat(clustering1, clustering2, greedy.colours = c("red"))
    # optimal bi-graph without greedy algorithm
    flatVSflat(clustering1, clustering2, greedy = FALSE)


    # simulated data; only one supercluster
    clustering1 <- c(rep(1, 5), rep(2, 10), rep(3, 10))
    clustering2 <- c(rep(1:4, 5), rep(1, 5))
    flatVSflat(clustering1, clustering2, greedy.colours = c("red"))

    
}
\keyword{ clustering comparison }

