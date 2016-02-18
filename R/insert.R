insert <- function (vect, position, value) {
    
    vect <- as.vector(vect)
    value <- as.vector(value)
    position <- as.numeric(position)
    if (length(position) > 1) {stop("Position must be a number...")}
    if (position == 0) {stop("Position must be a positive integer")}
    if (position == 1) {new.vector <- c(value, vect)}
    else {
        if (position > length(vect)) {
            new.vector <- vect
            new.vector[position:(position + length(value) - 1)] <- value
        }  else {
            new.vector <- c(vect[1:(position - 1)], value,
            vect[position:length(vect)])
        } # end ELSE
    } # end ELSE
    
    return(new.vector)
}
