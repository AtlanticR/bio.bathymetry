
bathymetry.parameters = function(DS="bio.bathymetry", p=NULL, resolution="canada.east.highres", nc=1 ) {

  if ( is.null(p) ) p=list()
  if ( !exists("project.name", p) ) p$project.name=DS

  if (DS=="bio.bathymetry"){
    p$project.root = project.datadirectory( p$project.name )
    p$libs = bioLibrary( "bio.base", "bio.utilities", "bio.bathymetry", "bio.coastline", "bio.polygons", "bio.conker", "conker" )
    p$libs = c( p$libs, RLibrary( c( "rgdal", "maps", "mapdata", "maptools", "lattice", "parallel", 
      "geosphere", "sp", "raster", "colorspace", "splancs" ) ) )
    # default (= only supported resolution of 0.5 km discretization)  .. do NOT change
    # use "complete" to project/downscale/upscale onto other grids/resolutions
    p = spatial_parameters( type=resolution, p=p )
    p = spatial_parameters(p)  # load defaults
    # cluster definition
    p$clusters = rep( "localhost", nc )
    return(p)
  }

  if (DS=="conker") {

    p$libs = RLibrary( c( p$libs, "conker" ) ) # required for parallel

    p$conker_rsquared_threshold = 0.3 # lower threshold
    p$conker_distance_prediction = 5 # this is a half window km
    p$conker_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    
    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  # fractions of median distance scale to try in local block search

    if (!exists("conker_local_modelengine", p)) p$conker_local_modelengine="inla"  # currently the perferred approach 

    if ( p$conker_local_modelengine =="gaussianprocess2Dt" ) {
      # too slow to use right now
      p$fields.cov.function = "stationary.cov"  #
      p$fields.cov.function = "stationary.taper.cov"  # Wendland tapering
      if (!exists("fields.Covariance", p)) p$fields.Covariance="Exponential" # note that "Rad.cov" is TPS
      if (!exists("fields.cov.args", p) & p$fields.Covariance=="Matern") {
        if (!exists("fields.nu", p)) p$fields.nu=0.5  # note: this is the smoothness or shape parameter (fix at 0.5 if not calculated or given -- exponential)   
        p$fields.cov.args=list( Covariance=p$fields.Covariance, smoothness=p$fields.nu ) # this is exponential covariance 
      }
    } else if (p$conker_local_modelengine == "kernel.density") {
      # ~ 3.25 days hr with 68, 3 Ghz cpus on beowulf using kernel.density method, bigmemory-filebacked jc: 2016 
      # ~ 20 hr with 8, 3.2 Ghz cpus on thoth using kernel.density method RAM based jc: 2016
      p$conker_local_modelformula = NULL # no need for formulae for kernel.density
    
    } else if (p$conker_local_modelengine == "gam") {
      # GAM are overly smooth .. adding more knots might be good but speed is the cost .. k=50 to 100 seems to work nicely
      ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
      p$conker_local_modelformula = formula( 
        z ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=100, bs="ts") )  
      p$conker_local_model_distanceweighted = TRUE  
      p$conker_local_family = bio.conker::log_gaussian_offset(1000)
    
    } else if ( p$conker_local_modelengine == "bayesx" ) {
    
      ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
      p$conker_local_modelformula = formula(z ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te") )  # more detail than "gs" .. "te" is preferred
      p$conker_local_model_bayesxmethod="MCMC"  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..
      p$conker_local_model_distanceweighted = TRUE  
      p$conker_local_family = bio.conker::log_gaussian_offset(1000)
      p$conker_local_family_bayesx ="gaussian"

    } else if (p$conker_local_modelengine == "inla" ){
      
      # old method .. took a month to finish .. results look good but very slow
      ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
      p$conker_local_family = bio.conker::log_gaussian_offset(1000)
      p$inla_family = "gaussian"
      p$inla.alpha = 0.5 # bessel function curviness .. ie "nu"
      p$conker_local_modelformula = formula( z ~ -1 + intercept + f( spatial.field, model=SPDE ) ) # SPDE is the spatial covaria0nce model .. defined in conker___inla
      p$conker.posterior.extract = function(s, rnm) {
        # rnm are the rownames that will contain info about the indices ..
        # optimally the grep search should only be done once but doing so would
        # make it difficult to implement in a simple structure/manner ...
        # the overhead is minimal relative to the speed of modelling and posterior sampling
        i_intercept = grep("intercept", rnm, fixed=TRUE ) # matching the model index "intercept" above .. etc
        i_spatial.field = grep("spatial.field", rnm, fixed=TRUE )
        return(  s$latent[i_intercept,1] + s$latent[ i_spatial.field,1] )
      }
    
    } else {

      message( "The specified conker_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }

    # other options might work depending upon data density but GP are esp slow .. too slow for bathymetry
    if (!exists("conker_variogram_method", p)) p$conker_variogram_method = "fast"
    
    p$variables = list( Y="z", LOCS=c("plon", "plat") ) 
    
    p$n.min = 30 # n.min/n.max changes with resolution
    p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time
  
    # if not in one go, then the value must be reconstructed from the correct elements:
    p$conker_sbox = conker_db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km

    p$conker_nonconvexhull_alpha = 50  # radius in distance units (km) to use for determining boundaries
    p$conker_theta = p$pres # FFT kernel bandwidth (SD of kernel) required for method "harmonic.1/kernel.density"

    p$conker_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$quantile_bounds = c(0.001, 0.999) # remove these extremes in interpolations

    return(p)
  }


}

