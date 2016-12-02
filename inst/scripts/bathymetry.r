
# Bathymetry data: processing bathymetry data with RINLA  .. no GMT dependency
# warning: this will take weeks as it is an iterative process

  p = bio.bathymetry::bathymetry.parameters()
  # p$clusters = c( rep( "nyx", nc ), rep ("tartarus", nc), rep("kaos", nc ) )


  ### -----------------------------------------------------------------
  # prepare data for modelling and prediction:: faster if you do this step on kaos (the fileserver)
  # also needs about 42 GB RAM, JC 2015
 if ( basedata.redo ) {
    # processing is  at "canada.east.highres" but temporarily drop to "canada.east to update raw data"
    bathymetry.db ( p=spatial_parameters( p=p, type="canada.east" ), DS="z.lonlat.rawdata.redo", additional.data=c("snowcrab", "groundfish") )
    bathymetry.db( p=p, DS="bathymetry.sthm.inputs.data.redo" )  # Warning: req ~ 15 min, 40 GB RAM (2015, Jae) data to model (with covariates if any)
    bathymetry.db( p=p, DS="bathymetry.sthm.inputs.prediction.redo" ) # i.e, pred locations (with covariates if any )
  }
  
  
  ### -----------------------------------------------------------------

  p$sthm_engine = "kernel.density"  # about 5 X faster than bayesx-mcmc method .. perferred for now
  # p$sthm_engine = "inla"  # about 5 X faster than bayesx-mcmc method .. perferred for now
  # p$sthm_engine = "gaussianprocess2Dt"  # too slow for the data density  
  # p$sthm_engine = "gam" # 2nd choice
  # p$sthm_eng ine = "bayesx" # too slow
  
  p = bio.bathymetry::bathymetry.parameters( p=p, DS="sthm" )

  landmask.redo= FALSE
  if (landmask.redo) {
    # re-run only if default resolution is altered ... very slow 1 hr?
    bathymetry.db( DS="landmasks.create", p=p ) 
  }

  # 2.2GB/process and 2.2 GB in parent 
  # boundary def takes too long .. too much data to process -- skip
  method = "bigmemory.ram"
  if (method=="bigmemory.ram") {
    # ~ 20 hr with 8, 3.2 Ghz cpus on thoth using kernel.density method jc: 2016
    p$clusters = rep("localhost", 8)
    data.call = 'bathymetry.db( p=p, DS="bathymetry.sthm.data" )'
    p = sthm( DATA=data.call, p=p, storage.backend="bigmemory.ram", boundary=FALSE )  
  }

  if (method=="bigmemory.filebacked") {
    # ~ 3.25 days hr with 68, 3 Ghz cpus on beowulf using kernel.density method jc: 2016
    p$clusters = c( rep( "nyx", 24 ), rep ("tartarus", 24), rep("kaos", 20 ) ) 
    p = sthm( DATA=data.call, p=p, storage.backend="bigmemory.filebacked", boundary=FALSE )  
  }


  # bring together stats and predictions and any other required computations: slope and curvature
  bathymetry.db( p=p, DS="bathymetry.sthm.finalize.redo" )
  # B = bathymetry( p=p, DS="bathymetry.sthm.finalize" )     # to see the assimilated data:

  ### -----------------------------------------------------------------
  # as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
  # if you want more, will need to add to the list and modify the selection criteria
  p$new.grids = c( "canada.east.highres", "canada.east", "SSE", "SSE.mpa" , "snowcrab")
  bathymetry.db( p=p, DS="complete.redo", grids.new=p$new.grids )
  bathymetry.db ( p=p, DS="baseline.redo" )   # filtering of areas and or depth to reduce file size, in planar coords only


  ### -----------------------------------------------------------------
  # "snowcrab" subsets do exist but are simple subsets of SSE
  # so only the lookuptable below is all that is important as far as bathymetry is concerned
  # both share the same initial domains + resolutions
  bathymetry.db( p=spatial_parameters( type="snowcrab" ), DS="lookuptable.sse.snowcrab.redo" ) # indices to map SSE to snowcrab


  ### -----------------------------------------------------------------
  # to recreate new polygons, run the following:
  bathyclines.redo = FALSE
  depths = c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 400, 500, 600, 700, 750, 800, 900,
               1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 )
  if( bathyclines.redo ) {
    # note these polygons are created at the resolution specified in p$spatial.domain ..
    # which by default is very high ("canada.east.highres" = 0.5 km .. p$pres ).
    # For lower one specify an appropriate p$spatial.domain
    plygn = isobath.db( p=p, DS="isobath.redo", depths=depths  )
  }



  ### -----------------------------------------------------------------
  plygn = isobath.db( p=p, DS="isobath", depths=depths  )

  coast = coastline.db( xlim=c(-68,-52), ylim=c(41,50), no.clip=TRUE )  # no.clip is an option for maptools::getRgshhsMap
  plot( coast, col="transparent", border="steelblue2" , xlim=c(-68,-52), ylim=c(41,50),  xaxs="i", yaxs="i", axes=TRUE )  # ie. coastline
  lines( plygn[ as.character(c( 100, 200, 300 ))], col="gray90" ) # for multiple polygons
  lines( plygn[ as.character(c( 500, 1000))], col="gray80" ) # for multiple polygons
  # plot( plygn, xlim=c(-68,-52), ylim=c(41,50))  # all isobaths commented as it is slow ..


  # or to get in projected (planar) coords as defined by p$spatial domain
  plygn = isobath.db( p=p, DS="isobath", depths=c(100) , crs=p$internal.crs ) # as SpatialLines
  plot(plygn)

  plygn_aslist = coordinates( plygn)
  plot( 0,0, type="n", xlim=c(-200,200), ylim=c(-200,200)  )
  lapply( plygn_aslist[[1]], points, pch="." )

  plygn_as_xypoints = coordinates( as( plygn, "SpatialPoints") )# ... etc...
  plot(plygn_as_xypoints, pch=".",  xaxs="i", yaxs="i", axes=TRUE)


  # a few plots :

  p = spatial_parameters( type="canada.east.highres", p=p )
  b = bathymetry.db( p=p, DS="complete" )

  vn = "z"
  u = b[[vn]] [ is.finite( b[[vn]]) ]
  if (length( u) > 0 ) out = x[u]

  mypalette = colorRampPalette(c("darkblue","blue3", "green", "yellow", "orange","red3", "darkred"), space = "Lab")(100)
  mybreaks = classIntervals( u, n=length(mypalette), style="quantile")$brks

  depths = c( 100, 200, 300, 400, 500 )
  plygn = isobath.db( p=p, DS="isobath", depths=depths  )

  sab = as( Polygon( coords=polygon.db( id="aoi.st.anns") ), "SpatialPolygons" )
  sab = spTransform( sab, crs( plygn) )

  sp.layout= list( sab, plygn )

  spplot( b, vn, col.regions=mypalette, main=vn, sp.layout=coastLayout, col="transparent" )


  #### New method -- "fast"(er) estimation of covariance function
  # using geoR .. most stable and flexible approach so far, uses ML methods
  # spBayes a little to unstable and slow



