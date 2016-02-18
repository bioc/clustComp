test_dyn.cross <- function() {
    checkEquals(dyn.cross(matrix(c(4,1,2,3),2,2)),2)
    checkEquals(dyn.cross(matrix(c(4,0,0,3),2,2)),0)
    checkEquals(dyn.cross(matrix(c(0,4,3,0),2,2)),12)
}