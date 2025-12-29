# This is a general utility function to find distance between two vectors of locations. Obtains a matrix. #

e2dist <-
function (loc1, loc2)
{
    i <- sort(rep(1:nrow(loc2), nrow(loc1)))
    dvec <- sqrt((loc1[, 1] - loc2[i, 1])^2 + (loc1[, 2] - loc2[i, 2])^2)
    matrix(dvec, nrow = nrow(loc1), ncol = nrow(loc2), byrow = F)
}
 