test_barycentre <- function() {
    checkEquals(barycentre(1:3),5/3)
    checkEquals(barycentre(1:3,1:3),7/3)
}
