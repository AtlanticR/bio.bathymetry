
# Bathymetry data: processing bathymetry data with RINLA  .. no GMT dependency
# warning: this will take weeks as it is an iterative process

  p = list( project.name = "bio.bathymetry" )

  p$project.root = project.datadirectory( p$project.name )

  p$libs = bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.coastline", "bio.polygons" )
  p$libs = c( p$libs, RLibrary( c( "rgdal", "maps", "mapdata", "maptools", "lattice", "parallel", "INLA", "gstat", "geoR",
    "geosphere", "sp", "raster", "colorspace" ,  "splancs", "fields",
    "bigmemory.sri", "synchronicity", "bigmemory", "biganalytics", "bigtabulate", "bigalgebra" ) ) )

  # default (= only supported resolution of 0.5 km discretization)  .. do NOT change
  # use "complete" to project/downscale/upscale onto other grids/resolutions
  p = spatial.parameters( type="canada.east.highres", p=p )

  p = spacetime.parameters(p)  # load defaults

  p$bathymetry.bigmemory.reset = FALSE

  # cluster definition
  nc = 1
  # nc = 5
  p$clusters = rep( "localhost", nc )
  # p$clusters = c( rep( "nyx", nc ), rep ("tartarus", nc), rep("kaos", nc ) )


  ### -----------------------------------------------------------------
  make.bathymetry.db = FALSE
  if (make.bathymetry.db) {
    # prepare data for modelling and prediction:: faster if you do this step on kaos (the fileserver)
    bathymetry.db ( p=spatial.parameters( type="canada.east", p=p ), DS="z.lonlat.rawdata.redo",
      additional.data=c("snowcrab", "groundfish") )
  }


  ### -----------------------------------------------------------------
  spatial.covariance.redo = FALSE
  if (spatial.covariance.redo) {
    p$clusters = c( rep( "nyx", 24 ), rep ("tartarus", 24), rep("kaos", 24 ) )
    # p$bathymetry.bigmemory.reset = TRUE   # reset needed if variables entering are changing (eg., addiing covariates with interpolation, etc)
    bathymetry.db( p=p, DS="covariance.spatial.redo" )
  }
  covSp = bathymetry.db( p=p, DS="covariance.spatial" )


  ### -----------------------------------------------------------------
  spatial.interpolation.redo = FALSE
  if (spatial.interpolation.redo) {
    # do not use all CPU's as INLA itself is partially run in parallel
    # RAM reqiurements are a function of data density and mesh density .. currently ~ 12 GB / run
    p$clusters = c( rep( "nyx", 5 ), rep ("tartarus", 5), rep("kaos", 5 ) )
    # p$bathymetry.bigmemory.reset = TRUE   # reset needed if variables entering are changing (eg., addiing covariates with interpolation, etc)
    # bathymetry.db( DS="landmasks.create", p=p ) # re-run only if default resolution is altered ... very slow 1 hr?
    bathymetry.db( p=p DS="spde.redo" )
    # to see the raw saved versions of the the results:
    # predSp = spacetime.db( p=p, DS="predictions.redo" )
    # statSp = spacetime.db( p=p, DS="statistics.redo" )
    # to see the assimilated data:
    # B = bathymetry( p=p, DS="bathymetry.spacetime.finalize" )
  }


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
  bathymetry.db( p=spatial.parameters( type="snowcrab" ), DS="lookuptable.sse.snowcrab.redo" ) # indices to map SSE to snowcrab


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

  p = spatial.parameters( type="canada.east.highres", p=p )
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



