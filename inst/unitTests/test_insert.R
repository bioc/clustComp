test_insert <- function() {
    checkEquals(insert(1:3,4,4), c(1:4))
    checkEquals(insert(1:3,4,4:6),c(1:6))
    checkEquals(length(insert(1:2,4,4)),4)
    checkTrue(any(is.na(insert(1:2,4,4))))
}
