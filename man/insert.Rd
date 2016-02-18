\name{insert}
\alias{insert}
\title{ Insert a set of values at a given position of a vector}
\description{
\code{insert} introduces a vector at a given position of another vector, 
displacing to the right all values from that position onwards. }
\usage{
insert(vect, position, value)
}
\arguments{
    \item{vect}{ the vector in which to insert additional values.} 
    \item{position}{ the position of vector \code{vect} at which to insert 
additional values.}
    \item{value}{ the values to be inserted in the vector \code{vect}.} 
}
\value{
    \item{new.vector}{ the vector resulting after inserting \code{value} at 
the position determined by \code{position}.}
}
\details{ The value of \code{position} does not need to be a number smaller 
than or equal to the length of \code{vect}, as the missing values will be 
denoted as NA. }
\author{
Aurora Torrente \email{aurora@ebi.ac.uk}
and
Alvis Brazma \email{brazma@ebi.ac.uk}
}
\examples{
    # simulated data
    v1 <- 1:10
    v2 <- insert(v1, 4, 0)
    v2 <- insert(v1, 1, rep(0, 5))
    v2 <- insert(v1, 11, "A")
    v2 <- insert(v1, 12, "A")
}
\keyword{ branch split }
