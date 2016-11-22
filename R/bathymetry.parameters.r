
bathymetry.parameters = function(DS="bio.bathymetry", p=NULL, resolution="canada.east.highres", nc=1 ) {

  if ( is.null(p) ) p=list()
  if ( !exists("project.name", p) ) p$project.name=DS

  if (DS=="bio.bathymetry"){
    p$project.root = project.datadirectory( p$project.name )
    p$libs = bioLibrary( "bio.base", "bio.utilities", "bio.bathymetry", "bio.coastline", "bio.polygons", "bio.spacetime" )
    p$libs = c( p$libs, RLibrary( c( "rgdal", "maps", "mapdata", "maptools", "lattice", "parallel", "INLA", "gstat", "geoR",
      "geosphere", "sp", "raster", "colorspace" ,  "splancs", "fields",  "ff", "ffbase" ) ) )
    # default (= only supported resolution of 0.5 km discretization)  .. do NOT change
    # use "complete" to project/downscale/upscale onto other grids/resolutions
    p = spacetime_parameters( type=resolution, p=p )
    p = spacetime_parameters(p)  # load defaults
    # cluster definition
    p$clusters = rep( "localhost", nc )
    return(p)
  }

  if (DS=="bio.bathymetry.spacetime") {

    p$spacetime_rsquared_threshold = 0.3 # lower threshold
    p$spacetime_distance_prediction = 5 # this is a half window km
    p$spacetime_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  # fractions of median distance scale (dist.max, dist.min)/2 to try in local block search


    if ( p$spacetime_engine =="gaussianprocess2Dt" ) {
      p$fields.cov.function = "stationary.cov"  #
      # p$fields.cov.function = "stationary.taper.cov"  # Wendland tapering
      if (!exists("fields.Covariance", p)) p$fields.Covariance="Exponential" # note that "Rad.cov" is TPS
      if (!exists("fields.cov.args", p) & p$fields.Covariance=="Matern") {
        if (!exists("fields.nu", p)) p$fields.nu=0.5  # note: this is the smoothness or shape parameter (fix at 0.5 if not calculated or given -- exponential)   
        p$fields.cov.args=list( Covariance=p$fields.Covariance, smoothness=p$fields.nu ) # this is exponential covariance 
      }
    }

    if (0) {
      # GAM are overly smooth
      p$spacetime_engine = "gam" # see model form in spacetime.r (method="xyts")
      p$spacetime_engine_modelformula = formula( 
        z ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=50, bs="ts") )  
      # p$spacetime_engine_modelformula = formula( z ~ -1 + intercept + f( spatial.field, model=SPDE ) ) # SPDE is the spatial covaria0nce model .. defined in spacetime___inla
      
      p$spacetime_engine = "bayesx"
      p$spacetime_engine_modelformula = formula(z ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te") )  # more detail than "gs" .. "te" is preferred
      p$bayesx.method="MCMC"  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..

    }

    p$spacetime_engine = "kernel.density"  # about 5 X faster than bayesx-mcmc method
    p$spacetime_engine_modelformula = NA # no need for formulae for kernel.density

    p$spacetime_model_distance_weighted = TRUE  # not required for kernel.density

    p$spacetime_covariate_modeltype="gam"
    p$spacetime_covariate_modelformula = p$spacetime_engine_modelformula

    p$variables = list( Y="z", LOCS=c("plon", "plat") ) 

    ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
    log_gaussian_offset = function(offset=0) {
      structure(list(
        linkfun = function(mu) log(mu + offset), 
        linkinv = function(eta) exp(eta) - offset,
        mu.eta = function(eta) NA, 
        valideta = function(eta) TRUE, 
        name = paste0("logexp(", offset, ")") ),
        class = "link-glm" )
    }
    p$spacetime_family = log_gaussian_offset(1000)
    
    p$dist.max = 50 # length scale (km) of local analysis .. for acceptance into the local analysis/model
    p$dist.min = 2 # lower than this .. subsampling occurs
    p$n.min = 30 # n.min/n.max changes with resolution: at p$pres=0.25, p$dist.max=25: the max count expected is 40000
    p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time
  
    # if not in one go, then the value must be reconstructed from the correct elements:
    p$sbbox = spacetime_db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km

    p$non_convex_hull_alpha = 50  # radius in distance units (km) to use for determining boundaries
    p$theta = p$pres # FFT kernel bandwidth (SD of kernel) required for method "harmonic.1/kernel.density"

    p$spacetime.noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$quantile_bounds = c(0.001, 0.999) # remove these extremes in interpolations

    return(p)
  }


}

