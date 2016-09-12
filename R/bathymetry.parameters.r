
bathymetry.parameters = function(DS="bio.bathymetry", p=NULL, resolution="canada.east.highres", nc=1 ) {

  if ( is.null(p) ) p=list()
  if ( !exists("project.name", p) ) p$project.name=DS

  if (DS=="bio.bathymetry"){
    p$project.root = project.datadirectory( p$project.name )
    p$libs = bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.coastline", "bio.polygons" )
    p$libs = c( p$libs, RLibrary( c( "rgdal", "maps", "mapdata", "maptools", "lattice", "parallel", "INLA", "gstat", "geoR",
      "geosphere", "sp", "raster", "colorspace" ,  "splancs", "fields",  "bigmemory" ) ) )
    # default (= only supported resolution of 0.5 km discretization)  .. do NOT change
    # use "complete" to project/downscale/upscale onto other grids/resolutions
    p = spatial.parameters( type=resolution, p=p )
    p = spacetime.parameters(p)  # load defaults
    # cluster definition
    p$clusters = rep( "localhost", nc )
    return(p)
  }

  if (DS=="bio.bathymetry.spacetime") {
    p$rootdir = file.path( p$project.root, "spacetime" )
    p$fn.P =  file.path( p$rootdir, paste( "spacetime", "predictions", p$spatial.domain, "rdata", sep=".") )
    p$fn.S =  file.path( p$rootdir, paste( "spacetime", "statistics", p$spatial.domain, "rdata", sep=".") )
    p$fn.results.covar =  file.path( p$rootdir, paste( "spatial", "covariance", p$spatial.domain, "rdata", sep=".") )
    p$variogram.engine = "gstat"  # "geoR" seg faults frequently ..
    p$dist.mwin = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$upsampling = c( 1.1, 1.2, 1.5, 2 )  # local block search fractions
    p$downsampling = c( 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2 ) # local block search fractions  -- need to adjust based upon data density
    p$mesh.boundary.resolution = 150
    p$mesh.boundary.convex = -0.025
    p$variables = list( Y="z", LOCS=c("plon", "plat") )
    p$spacetime.link = function( X ) { log(X + 1000) }  ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one -order of magnitude
    p$spacetime.invlink = function( X ) { exp(X) - 1000 }
    p$dist.max = 100 # length scale (km) of local analysis .. for acceptance into the local analysis/model
    p$dist.min = 75 # lower than this .. subsampling occurs
    p$dist.pred = 0.95 # % of dist.max where **predictions** are retained (to remove edge effects)
    p$n.min = 30 # n.min/n.max changes with resolution: at p$pres=0.25, p$dist.max=25: the max count expected is 40000
    p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time
    p$expected.range = 50 #+units=km km , with dependent var on log scale
    p$expected.sigma = 1e-1  # spatial standard deviation (partial sill) .. on log scale
    p$spatial.field.name = "spatial.field"  # name used in formula to index the spatal random field
    p$modelformula = formula( z ~ -1 + intercept + f( spatial.field, model=SPDE ) ) # SPDE is the spatial covariance model .. defined in spacetime.interpolate.inla.local (below)
    p$spacetime.family = "gaussian"
    p$spacetime.outputs = c( "predictions.projected", "statistics" ) # "random.field", etc.
    p$statsvars = c("range", "range.sd", "spatial.error", "observation.error")
    # if not in one go, then the value must be reconstructed from the correct elements:
    p$sbbox = spacetime.db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km
    p$nPreds = p$nplons * p$nplats
     p$spacetime.posterior.extract = function(s, rnm) {
      # rnm are the rownames that will contain info about the indices ..
      # optimally the grep search should only be done once but doing so would
      # make it difficult to implement in a simple structure/manner ...
      # the overhead is minimal relative to the speed of modelling and posterior sampling
      i_intercept = grep("intercept", rnm, fixed=TRUE ) # matching the model index "intercept" above .. etc
      i_spatial.field = grep("spatial.field", rnm, fixed=TRUE )
      return(  s$latent[i_intercept,1] + s$latent[ i_spatial.field,1] )
    }
    return(p)
  }


}

