
bathymetry.figures = function( DS=NULL, p=NULL, zrange=NULL ) {

  bioLibrary( "bio.coastline" )

  if ( DS=="predictions" ) {

    p = spacetime.db( p=p, DS="bigmemory.filenames" )
    P = attach.big.matrix(p$descriptorfile.P , path=p$tmp.datadir )
    pps  =  expand.grid( plons=p$plons, plats=p$plats)
    p$spatial.domain="canada.east"  # force isobaths to work in levelplot
    if (is.null( zrange)) {
      datarange = log( c( 5, 4000 ))
    } else {
      datarange = log( zrange)
    }

    dr = seq( datarange[1], datarange[2], length.out=100)
    oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions" )
    # oc= 1:nrow(P[])
    print(
      levelplot( log( P[oc,2] ) ~ plons + plats, pps[oc,], aspect="iso", main=NULL, at=dr, col.regions=rev(color.code( "seis", dr)) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE),
        panel = function(x, y, subscripts, ...) {
          panel.levelplot (x, y, subscripts, aspect="iso", rez=c(1,1), ...)
          sp.lines( isobath.db( p=p, DS="isobath", depths=c(200, 400 ), crs=p$internal.crs ), col = "gray80", cex=0.1 )
          sp.lines( coastline.db( p=p, crs=p$internal.crs ), col = "steelblue", cex=0.1 )
        }
    ) )
  }

  if ( DS=="predictions.error" ) {
    p = spacetime.db( p=p, DS="bigmemory.filenames" )
    P = attach.big.matrix(p$descriptorfile.P , path=p$tmp.datadir )
    pps  =  expand.grid( plons=p$plons, plats=p$plats)
    p$spatial.domain="canada.east"  # force isobaths to work in levelplot
    if (is.null( zrange)) {
      datarange = log( c( 2, 50 ))
    } else {
      datarange = log( zrange)
    }

    dr = seq( datarange[1], datarange[2], length.out=100)
    oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions" )
    print(
     levelplot( log( P[oc,3] ) ~ plons + plats, pps[oc,], aspect="iso", main=NULL, at=dr, col.regions=rev(color.code( "seis", dr)) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE),
        panel = function(x, y, subscripts, ...) {
          panel.levelplot (x, y, subscripts, aspect="iso", rez=c(1,1), ...)
          sp.lines( isobath.db( p=p, DS="isobath", depths=c(200, 400 ), crs=p$internal.crs ), col = "gray80", cex=0.1  )
          sp.lines( coastline.db( p=p, crs=p$internal.crs ), col = "steelblue", cex=0.1 )
        }
     ))
  }


  if ( DS=="statistics" ) {
    p = spacetime.db( p=p, DS="bigmemory.filenames" )
    S = attach.big.matrix(p$descriptorfile.S , path=p$tmp.datadir )
    Sloc = attach.big.matrix(p$descriptorfile.Sloc , path=p$tmp.datadir )
    p$spatial.domain="canada.east"  # force isobaths to work in levelplot
    if (is.null( zrange)) {
      datarange = log( c( 5, 800 ))
    } else {
      datarange = log( zrange)
    }


    dr = seq( datarange[1], datarange[2], length.out=150)
    oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="statistics" )
    print( levelplot( log(S[oc,1])  ~ Sloc[oc,1] + Sloc[oc,2] , aspect="iso", at=dr, col.regions=color.code( "seis", dr) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE), cex=2,
      panel = function(x, y, subscripts, ...) {
        panel.levelplot (x, y, subscripts, aspect="iso", rez=c(5,5), ...)
        sp.lines( isobath.db( p=p, DS="isobath", depths=c( 200, 400 ), crs=p$internal.crs  ), col = "gray80", cex=0.1 )
        sp.lines( coastline.db( p=p, crs=p$internal.crs ), col = "steelblue", cex=0.1 )
      }
   ) )
  }




}


