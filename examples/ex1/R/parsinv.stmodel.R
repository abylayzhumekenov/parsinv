#' Define a spatio-temporal model with specified type and manifolds
#' 
#' The function creates a vector of integers, containing the the code / id
#' for model type, temporal and spatial manifolds. The allowed types are
#' `("102", "121", "202", "220")`. Temporal manifolds must be one of `("R1", "S1")`,
#' which are coded as `(10, 11)`. For spatial manifolds, these are `("R2", "S2")`
#' and `(20, 21)`.
#' 
#' @param type model type
#' @param manifold.t temporal manifold
#' @param manifold.s spatial manifold
#' @return integer vector of [type, manifold.t, manifold.s]
#' @export
parsinv.stmodel.define = function(type, manifold.t, manifold.s){
    stopifnot(type %in% c("102", "121", "202", "220"),
              manifold.t %in% c("R1", "S1"),
              manifold.s %in% c("R2", "S2"))
    manifold.code = list(R1 = 10, S1 = 11, R2 = 20, S2 = 21)
    return(as.integer(c(type, manifold.code[[manifold.t]], manifold.code[[manifold.s]])))
}


#' Write the model info to a file
#' 
#' The function writes a model information into a binary file.
#' 
#' @param stmodel model
#' @param filename file to write
#' @export 
parsinv.stmodel.write = function(stmodel, filename){
    con = file(filename, "wb")
    writeBin(object = stmodel, con = con, size = 4)
    close(con = con)
}


#' Read the model from a file
#' 
#' The function reads a model information from a binary file.
#' 
#' @param filename file to read from
#' @return integer vector of [type, manifold.t, manifold.s]
#' @export
parsinv.stmodel.read = function(filename){
    con = file(filename, "rb")
    stmodel = readBin(con = con, what = "integer", n = 3)
    close(con = con)
    return(stmodel)
}
