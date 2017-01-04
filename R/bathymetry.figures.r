
bathymetry.figures = function( p=NULL, vn="z", datarange=NULL, logyvar=FALSE, isodepths = c( 100, 300, 500 ) ) {

  bioLibrary( "bio.coastline" )

  b = bathymetry.db( p=p, DS="complete" )
  if (logyvar) {
    b = b[ b[,vn]>0, ]
    b[,vn] = log(b[,vn])
  }

  oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions",internal.projection=p$internal.projection )
  if (length(oc) > 0) {
    b = b[oc,]
  }

  if (is.null( datarange)) datarange = quantile( b[,vn], probs=c(0.05, 0.95), na.rm=TRUE )
  dr = seq( datarange[1], datarange[2], length.out=100)

  print( 
    levelplot( b[,vn] ~ plon + plat, data=b[oc,], aspect="iso", main=NULL, 
      at=dr, col.regions=rev(color.code( "seis", dr)) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE),
      panel = function(x, y, subscripts, ...) {
        panel.levelplot (x, y, subscripts, aspect="iso", rez=c(1,1), ...)
        sp.lines( isobath.db( p=p, DS="isobath", depths=isodepths, crs=p$internal.crs ), col = "gray80", cex=0.1 )
        sp.lines( coastline.db( p=p, crs=p$internal.crs ), col = "steelblue", cex=0.1 )
      }
  ) )

}


