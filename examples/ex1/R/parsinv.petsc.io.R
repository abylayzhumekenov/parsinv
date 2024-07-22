#' Write a sparse matrix in binary CSR format with PETSc header
#'
#' The function saves a matrix as a binary file in compressed sparse row format.
#' The header also includes classid of the PETSc Mat object.
#' Ignore first 4 bytes, if you read the matrix from elsewhere.
#' Additionally, endian = "swap" for Linux machines.
#'
#' @param A sparse matrix
#' @param filename binary file to write
#' @export
parsinv.mat.write = function(A, filename){
    # encode the matrix
    A = as(A, "RsparseMatrix")
    A = as(A, "generalMatrix")
    x = list(classid = 1211216,
             nrows = nrow(A),
             ncols = ncol(A),
             nnz = length(A@x),
             nnz_row = diff(A@p),
             nnz_i = A@j,
             nnz_val = A@x)

    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$ncols), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz_row), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz_i), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.double(x$nnz_val), con = fwrite, size = 8, endian = "swap")
    close(con = fwrite)
}


#' Write a dense matrix in binary CSR format with PETSc header
#'
#' The function saves a matrix as a binary file in compressed sparse row format.
#' The header also includes classid of the PETSc Mat object.
#' Ignore first 4 bytes, if you read the matrix from elsewhere.
#' Additionally, endian = "swap" for Linux machines.
#'
#' @param Q dense matrix
#' @param filename binary file to write
#' @export
parsinv.dense.write = function(A, filename){
    # encode the matrix
    x = list(classid = 1211216,
             nrows = nrow(A),
             ncols = ncol(A),
             nnz = length(A),
             nnz_row = rep(ncol(A), nrow(A)),
             nnz_i = rep(1:ncol(A)-1, times = nrow(A)),
             nnz_val = c(t(A)))

    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$ncols), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz_row), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz_i), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.double(x$nnz_val), con = fwrite, size = 8, endian = "swap")
    close(con = fwrite)
}


#' Write a vector of doubles to a binary file with PETSc header
#'
#' The function saves a vector of doubles as a binary file.
#' The header also includes classid of the PETSc Vec object.
#' Ignore first 4 bytes, if you read the vector from elsewhere.
#' Additionally, endian = "swap" for Linux machines.
#'
#' @param y vector
#' @param filename binary file to write
#' @export
parsinv.vec.write = function(y, filename){
    # encode the vector
    x = list(classid = 1211214,
             nrows = length(y),
             val = drop(y))

    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.double(x$val), con = fwrite, size = 8, endian = "swap")
    close(con = fwrite)
}


#' Write a vector of integers to a binary file with PETSc header
#'
#' The function saves a vector of integers as a binary file.
#' The header also includes classid of the PETSc IS object.
#' Ignore first 4 bytes, if you read the vector from elsewhere.
#' Additionally, endian = "swap" for Linux machines.
#'
#' @param is index set (vector of integers)
#' @param filename binary file to write
#' @export
parsinv.is.write = function(is, filename){
    # encode the index set
    x = list(classid = 1211218,
             nrows = length(is),
             val = drop(is)-1)

    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$val), con = fwrite, size = 4, endian = "swap")
    close(con = fwrite)
}


#' Read a vector of doubles from a binary file with PETSc header
#'
#' The function loads a vector of doubles from a binary file.
#' The header also includes classid of the PETSc Vec object.
#' We ignore first 4 bytes, and "swap" the endianness for Linux machines.
#'
#' @param y vector
#' @param filename binary file to write
#' @return numeric vector
#' @export
parsinv.vec.read = function(filename){
    # decode the vector
    fread = file(filename, "rb")
    classid = readBin(con = fread, what = "integer", n = 1, endian = "swap")
    nrows = readBin(con = fread, what = "integer", n = 1, endian = "swap")
    val = readBin(con = fread, what = "double", n = nrows, endian = "swap")
    close(con = fread)
    return(val)
}


#' Write a named list of sparse matrices into binary files.
#'
#' The function saves sparse matrices stores in a named list into corresponding
#' binary files in a given directory.
#'
#' @param matrices a list of sparse matrices
#' @param dirname directory to save file in
#' @export
parsinv.mats.write = function(matrices, dirname){
    filenames = names(matrices)
    stopifnot(length(unique(filenames)) == length(matrices))
    for(i in seq_along(matrices)){
        parsinv.mat.write(matrices[[i]], paste0(dirname, "/", filenames[i]))
    }
}
