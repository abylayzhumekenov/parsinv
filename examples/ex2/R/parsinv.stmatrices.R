#' Create J matrices from the temporal mesh
#' 
#' The function generates temporal precision matrices for the Kronecker product.
#' 
#' @param stmodel model
#' @param tmesh temporal mesh (for now, only regular time points)
#' @return list of sparse J matrices
#' @export
parsinv.jmatrices = function(stmodel, tmesh){
    # compute the order
    alpha_t = stmodel[1] %/% 100 %% 10
    order = alpha_t
    
    # create a fem object
    tfem = fmesher::fm_fem(tmesh, order = order)
    n = nrow(tfem$c0)
    h = mean(diff(tmesh$loc))
    
    # 1st order in time
    if(stmodel[1] %in% c(102, 121)){
        J0 = tfem$c0
        J1 = Matrix::sparseMatrix(i = c(1, n), j = c(1, n), x = 0)
        J2 = tfem$g1
        
        # correction if not cyclic
        if(!tmesh$cyclic){
            J1[1,1] = J1[n,n] = 0.5
        }

        return(list(J0 = J0, J1 = J1, J2 = J2))
    }
    
    # 2nd order in time
    if(stmodel[1] %in% c(202, 220)){
        J0 = tfem$c0
        J1 = Matrix::sparseMatrix(i = c(1,1,2,2,n-1,n-1,n,n),
                                  j = c(1,2,1,2,n-1,n,n-1,n),
                                  x = 0)
        J2 = tfem$g1 * 2
        J3 = Matrix::sparseMatrix(i = c(1,1,2,2,n-1,n-1,n,n),
                                  j = c(1,2,1,2,n-1,n,n-1,n),
                                  x = 0)
        J4 = tfem$g2
        
        # corrections if not cyclic
        if(!tmesh$cyclic){
            J0[2,2] = tfem$c0[2,2] / 2
            J0[n-1,n-1] = tfem$c0[n-1,n-1] / 2
            
            J1[1:2,1:2] = J1[(n-1):n,(n-1):n] = c(5,-1,-1,5)/4
            
            J2[1,2] = J2[2,1] = tfem$g1[1,2]
            J2[2,2] = tfem$g1[2,2]
            J2[n-1,n] = J2[n,n-1] = tfem$g1[n-1,n]
            J2[n-1,n-1] = tfem$g1[n-1,n-1]
            
            J3[1:2,1:2] = J3[(n-1):n,(n-1):n] = c(2,-2,-2,2)/h^2
            
            J4[1,1] = tfem$g2[1,1] / 3
            J4[1,2] = J4[2,1] = tfem$g2[1,2] / 2
            J4[2,2] = tfem$g2[2,2] * 5/7
            J4[n,n] = tfem$g2[n,n] / 3
            J4[n-1,n] = J4[n,n-1] = tfem$g2[n-1,n] / 2
            J4[n-1,n-1] = tfem$g2[n-1,n-1] * 5/7
        }

        return(list(J0 = J0, J1 = J1, J2 = J2, J3 = J3, J4 = J4))
    }
}


#' Create G matrices from the spatial mesh
#' 
#' The function generates spatial matrices for the Kronecker product.
#' 
#' @param stmodel model
#' @param smesh spatial mesh
#' @return list of sparse G matrices 
#' @export
parsinv.gmatrices = function(stmodel, smesh){
    # compute the order
    alpha_t = stmodel[1] %/% 100 %% 10
    alpha_s = stmodel[1] %/% 10 %% 10
    alpha_e = stmodel[1] %/% 1 %% 10
    order = alpha_s * alpha_t + alpha_e
    
    # create a fem object
    sfem = fmesher::fm_fem(smesh, order = order)
    
    # save to another list
    gmatrices = list(G0 = sfem$c0)
    for(i in 1:order){
        gmatrices[[paste0("G",i)]] = sfem[[paste0("g",i)]]
    }
    
    return(gmatrices)
}


#' Create A projection matrices from temporal and spatial meshes
#' 
#' The function generate sparse projection matrices from temporal and spatial
#' mesh objects and observation location points. The projection matrices
#' are transposed (for partitioning purposes).
#' 
#' @param tmesh temporal mesh 
#' @param smesh spatial mesh 
#' @param tloc temporal locations 
#' @param sloc spatial locations
#' @return list of temporal and spatial projection matrices
#' @export
parsinv.amatrices = function(tmesh, smesh, tloc = NULL, sloc = NULL){
    if(is.null(tloc)){
        At = t(inla.spde.make.A(tmesh))
    } else {
        At = t(inla.spde.make.A(tmesh, loc = tloc))
    }
    if(is.null(sloc)){
        As = t(inla.spde.make.A(smesh))
    } else {
        As = t(inla.spde.make.A(smesh, loc = sloc))
    }
    return(list(At = At, As = As))
}
