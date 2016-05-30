block.data.frame = function( X, byvars, function.block=block.mean ) {
  
    vnames0 = names( X ) 
    vnames = setdiff( vnames0, byvars )
    nv = length( vnames )

    bk = apply( X[,byvars], 1, paste0, collapse="_" )
    bku = unique( bk )
    nbku = length( bku )
    bkq = matrix( unlist( strsplit( bk, "_") ), nrow=nbku, byrow=TRUE )
    # bkq = 
    out = matrix(NA, nrow=nbku, ncol=nv )  

    for ( i in 1:nbku ) {
      rn = which( bk == bku[i] )
      for (j in 1:nv) {
        out[i,j] = function.block( X[rn,j] )  
      }
    }
    out = as.data.frame( cbind( bkq, out ) ) 
    
    out[,xi] = as.numeric(as.character( out[,xi] ))
     out[,yi] = as.numeric(as.character( out[,yi] ))
     out = out[ which( is.finite( out[, zi] )) ,]  
     names(out) = c(byvars, vars )
    return( out )
 }

