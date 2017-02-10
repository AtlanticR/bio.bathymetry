
# Bathymetry 
# warning: this will take weeks as it is an iterative process


if ( basedata.redo ) {
  p = bio.bathymetry::bathymetry.parameters() 
  bathymetry.db( p=p, DS="z.lonlat.rawdata.redo" ) # needs about 42 GB RAM, JC 2015
  bathymetry.db( p=p, DS="lbm.inputs.redo" )  # Warning: req ~ 15 min, 40 GB RAM (2015, Jae) data to model (with covariates if any)
}
  

### -----------------------------------------------------------------
# Spatial interpolation using lbm
# Total at "superhighres": 30GB --~3.5 GB/process and 4 GB in parent for fft; gam method requires more ~ 2X
# boundary def takes too long .. too much data to process -- skip
# "highres": ~ 20 hr with 8, 3.2 Ghz cpus on thoth using fft method jc: 2016 or ~ 6 hr on hyperion
# "superhighres": ~ 40hr with 8 cpu on thoth for "fft"
# "superhighres": ~ 61hr with 8 cpu on thoth for "krige" method  
#   -- looks to be the best in performance/quality; req ~5 GB per process req

p = bio.bathymetry::bathymetry.parameters() # reset to defaults
p$lbm_local_modelengine = "krige"  
p$storage.backend="bigmemory.ram"
p = bio.bathymetry::bathymetry.parameters( p=p, DS="lbm" )
# p$clusters = rep("localhost",  detectCores() )

DATA='bathymetry.db( p=p, DS="lbm.inputs" )'
lbm( p=p, tasks=c( "initiate" ), DATA=DATA ) # a few minutes
lbm( p=p, tasks=c( "stage1" ) )  #60hrs
lbm( p=p, tasks=c( "stage2" ) )  #10hrs ?
lbm( p=p, tasks=c( "stage3" ) )  #1hrs ?
lbm( p=p, tasks=c( "save" ) )

  # to view progress in terminal:
  # watch -n 120 cat /home/jae/bio.data/bio.bathymetry/modelled/t/canada.east.superhighres/lbm_current_status

  # to view maps from an external R session:
  # lbm(p=p, tasks="debug_pred_static_map", vindex=1)
  # lbm(p=p, tasks="debug_pred_static_log_map", vindex=1)
  # lbm(p=p, tasks="debug_pred_dynamic_map", vindex=1)
  # lbm(p=p, tasks="debug_stats_map", vindex=1)


# bring together stats and predictions and any other required computations: slope and curvature
# and then regrid/warp as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
# .. still uses about 30-40 GB as the base layer is "superhighres" .. 
# if parallelizing .. use different servers than local nodes
bathymetry.db( p=p, DS="complete.redo" ) # finalise at diff resolutions 15 min ..
bathymetry.db( p=p, DS="baseline.redo" )  # coords of areas of interest ..filtering of areas and or depth to reduce file size, in planar coords only



### -----------------------------------------------------------------
# to recreate new polygons, run the following:
bathyclines.redo = FALSE
depthsall = c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 750, 800, 900,
             1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 )
if( bathyclines.redo ) {
  # note these polygons are created at the resolution specified in p$spatial.domain ..
  # which by default is very high ("canada.east.highres" = 0.5 km .. p$pres ).
  # For lower one specify an appropriate p$spatial.domain
  options(max.contour.segments=10000) # required if superhighres is being used
  for (g in c("canada.east.superhighres", "canada.east.highres", "canada.east", "SSE", "SSE.mpa", "snowcrab")) {
    print(g)
    p = bio.bathymetry::bathymetry.parameters(resolution=g) 
    if( g=="snowcrab") depths = c( 10, 20, 50, 75, 100, 200, 250, 300, 350 )  # by definition .. in filter.bathymetry
    if( g=="SSE") depths = depthsall[ depthsall < 801] # by definition
    if( g=="SSE.mpa") depths = depthsall[depthsall<2001]  # by definition
    if( grepl( "canada.east", g)) depths = depthsall  
    plygn = isobath.db( p=p, DS="isobath.redo", depths=depths  )
  }
}


### -----------------------------------------------------------------
# some test plots

p = bio.bathymetry::bathymetry.parameters(resolution="canada.east") # reset to lower resolution 
  
plygn = isobath.db( p=p, DS="isobath", depths=depths  )

coast = coastline.db( xlim=c(-75,-52), ylim=c(41,50), no.clip=TRUE )  # no.clip is an option for maptools::getRgshhsMap
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

p = bio.bathymetry::bathymetry.parameters(resolution="canada.east.highres")
bathymetry.figures( p=p, vn="z", logyvar=TRUE )


p = bio.bathymetry::bathymetry.parameters(resolution="canada.east.superhighres")
bathymetry.figures( p=p, vn="z", logyvar=TRUE )


p = bio.bathymetry::bathymetry.parameters(resolution="snowcrab")
bathymetry.figures( p=p, vn="z", logyvar=TRUE )

