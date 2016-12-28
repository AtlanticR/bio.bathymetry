
bathymetry.parameters = function(DS="bio.bathymetry", p=NULL, resolution=NULL ) {

  if ( is.null(p) ) p=list()
  if ( !exists("project.name", p) ) p$project.name=DS

  if (DS=="bio.bathymetry"){
    p$project.root = project.datadirectory( p$project.name )
    p$libs = bioLibrary( "bio.base", "bio.utilities", "bio.bathymetry", "bio.coastline", "bio.polygons", "bio.spacetime", "hivemod" )
    p$libs = c( p$libs, RLibrary( c( "rgdal", "maps", "mapdata", "maptools", "lattice", "parallel", 
      "geosphere", "sp", "raster", "colorspace", "splancs" ) ) )
 
    p$default.spatial.resolution = "canada.east.superhighres"
    if (is.null( resolution )) resolution=p$default.spatial.resolution
    p = spatial_parameters(p=p, type=resolution )  # default (= only supported resolution of 0.2 km discretization)  .. do NOT change
    return(p)
  }


  if (DS=="hivemod") {

    p$libs = RLibrary( c( p$libs, "hivemod" ) ) # required for parallel
    p$clusters = rep("localhost", detectCores() )
    p$storage.backend="bigmemory.ram"
    
    p$boundary = FALSE
    p$depth.filter = FALSE # need data above sea level to get coastline
    p$hivemod_nonconvexhull_alpha = 20  # radius in distance units (km) to use for determining boundaries
    p$hivemod_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$hivemod_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
    
    p$hivemod_rsquared_threshold = 0.1 # lower threshold
    p$hivemod_distance_prediction = 7.5 # this is a half window km
    p$hivemod_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$hivemod_distance_scale = 25 # km ... approx guess of 95% AC range 
    p$hivemod_distance_min = 1
    p$hivemod_distance_max = 50 # never go beyond the min and max range ( cpu/ram and time consideration are part of it but mostly what is physically reasonable)

    
    p$n.min = 30 # n.min/n.max changes with resolution
    p$n.max = 5000 # numerical time/memory constraint -- anything larger takes too much time
    # other options might work depending upon data density but GP are esp slow .. too slow for bathymetry
    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  # fractions of median distance scale to try in local block search
 
    p$variables = list( Y="z", LOCS=c("plon", "plat") ) 

    if (!exists("hivemod_variogram_method", p)) p$hivemod_variogram_method = "fast"
   
    if (!exists("hivemod_local_modelengine", p)) p$hivemod_local_modelengine="krige"  # currently the perferred approach ; fft is faster but smoothing is too erratic

    if ( p$hivemod_local_modelengine %in% c("krige" )) { 
      p$hivemod_krige_engine="fields" # faster than gstat
      p$hivemod_local_family = hivemod::log_gaussian_offset(1000)

    } else if ( p$hivemod_local_modelengine =="gaussianprocess2Dt" ) {
      # too slow to use right now

      p$fields.cov.function = "stationary.cov"  #
      p$fields.cov.function = "stationary.taper.cov"  # Wendland tapering
      if (!exists("fields.Covariance", p)) p$fields.Covariance="Exponential" # note that "Rad.cov" is TPS
      if (!exists("fields.cov.args", p) & p$fields.Covariance=="Matern") {
        if (!exists("fields.nu", p)) p$fields.nu=0.5  # note: this is the smoothness or shape parameter (fix at 0.5 if not calculated or given -- exponential)   
        p$fields.cov.args=list( Covariance=p$fields.Covariance, smoothness=p$fields.nu ) # this is exponential covariance 
      }
      p$hivemod_local_family = hivemod::log_gaussian_offset(1000)

    } else if (p$hivemod_local_modelengine == "fft") {
      # ~ 3.25 days hr with 68, 3 Ghz cpus on beowulf using fft method, bigmemory-filebacked jc: 2016 
      # ~ 14 hrs with 8, 3.2 Ghz cpus on thoth; 1 GB per process and a total of 6 GB usage;  method RAM based jc: 2016
      # 12 hrs to complete stage 1 on hyperion 
      # ~ 5.5 hr on hyperion 
      # definitely a cleaner (not overly smoothed) image than a GAM
      # NOTE that  p$hivemod_lowpass_phi and  p$hivemod_lowpass_nu are very critical choices
      p$hivemod_local_family = hivemod::log_gaussian_offset(1000)

      p$hivemod_fft_filter = "lowpass" # only act as a low pass filter .. depth has enough data for this. Otherwise, use: 
      # p$hivemod_fft_filter = "spatial.process" to ~ krige
      p$hivemod_lowpass_phi = 0.5 # low pass FFT filter range .. 0.5 seems to be optimal (by visual inspection)
      p$hivemod_lowpass_nu = 0.5 # this is exponential covar

    } else if (p$hivemod_local_modelengine == "twostep") {

      stop( "Makes no sense to use two step as there is no time series")
    
    } else if (p$hivemod_local_modelengine == "gam") {
      # GAM are overly smooth .. adding more knots might be good but speed is the cost .. k=50 to 100 seems to work nicely
      # timings: 
      # 14 hrs on hyperion with 100 knots
      ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
      p$hivemod_local_modelformula = formula( 
        z ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=100, bs="ts") )  
      p$hivemod_local_model_distanceweighted = TRUE  
      p$hivemod_local_family = hivemod::log_gaussian_offset(1000)
    
    } else if ( p$hivemod_local_modelengine == "bayesx" ) {
    
      ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
      p$hivemod_local_modelformula = formula(z ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te") )  # more detail than "gs" .. "te" is preferred
      p$hivemod_local_model_bayesxmethod="MCMC"  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..
      p$hivemod_local_model_distanceweighted = TRUE  
      p$hivemod_local_family = hivemod::log_gaussian_offset(1000)
      p$hivemod_local_family_bayesx ="gaussian"

    } else if (p$hivemod_local_modelengine == "inla" ){
      
      # old method .. took a month to finish .. results look good but very slow
      ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
      p$hivemod_local_family = hivemod::log_gaussian_offset(1000)
      p$inla_family = "gaussian"
      p$inla.alpha = 0.5 # bessel function curviness .. ie "nu"
      p$hivemod_local_modelformula = formula( z ~ -1 + intercept + f( spatial.field, model=SPDE ) ) # SPDE is the spatial covaria0nce model .. defined in hivemod___inla
      p$hivemod.posterior.extract = function(s, rnm) {
        # rnm are the rownames that will contain info about the indices ..
        # optimally the grep search should only be done once but doing so would
        # make it difficult to implement in a simple structure/manner ...
        # the overhead is minimal relative to the speed of modelling and posterior sampling
        i_intercept = grep("intercept", rnm, fixed=TRUE ) # matching the model index "intercept" above .. etc
        i_spatial.field = grep("spatial.field", rnm, fixed=TRUE )
        return(  s$latent[i_intercept,1] + s$latent[ i_spatial.field,1] )
      }
    
    } else {

      message( "The specified hivemod_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }


    return(p)
  }


}

