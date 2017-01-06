
# Bathymetry 
# warning: this will take weeks as it is an iterative process


if ( basedata.redo ) {
  p = bio.bathymetry::bathymetry.parameters() 
  bathymetry.db( p=p, DS="z.lonlat.rawdata.redo" ) # needs about 42 GB RAM, JC 2015
  bathymetry.db( p=p, DS="landmasks.create" ) # re-run only if default resolution .. v. slow ... currently using sp::over ... replace with point.in.polygon (TODO)
  bathymetry.db( p=p, DS="lbm.inputs.redo" )  # Warning: req ~ 15 min, 40 GB RAM (2015, Jae) data to model (with covariates if any)
}
  

### -----------------------------------------------------------------
# Spatial interpolation using lbm
# Total at "superhighres": 30GB --~3.5 GB/process and 4 GB in parent for fft; gam method requires more ~ 2X
# boundary def takes too long .. too much data to process -- skip
# "highres": ~ 20 hr with 8, 3.2 Ghz cpus on thoth using fft method jc: 2016 or ~ 6 hr on hyperion
# "superhighres": ~ 40hr with 8 cpu on thoth
# "superhighres": ~ 60hr with 8 cpu on thoth for "krige" method  
#   -- looks to be the best in performance/quality; req ~5 GB per process req

p = bio.bathymetry::bathymetry.parameters() # reset to defaults
p$lbm_local_modelengine = "krige"  
p$storage.backend="bigmemory.ram"
p = bio.bathymetry::bathymetry.parameters( p=p, DS="lbm" )
# p$clusters = rep("localhost",  detectCores() )

DATA='bathymetry.db( p=p, DS="lbm.inputs" )'
p = lbm( p=p, tasks=c( "initiate" ), DATA=DATA )
p = lbm( p=p, tasks=c( "stage1", "stage2", "stage3" ) )
p = lbm( p=p, tasks=c( "save" ) )

# bring together stats and predictions and any other required computations: slope and curvature
bathymetry.db( p=p, DS="lbm.finalize.redo" ) 
# B = bathymetry( p=p, DS="lbm.finalize" )     # to see the assimilated data:


# as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
# if you want more, will need to add to the list and modify the selection criteria
# .. still uses about 30-40 GB as the base layer is "superhighres" .. 
# if parallelizing .. use different servers than local nodes
  # this requires "raster" (it is possible to use fields and be a bit faster but this is simpler for now .. do  not forget to take care of the different projections)
p$new.grids = c( "canada.east.superhighres", "canada.east.highres", "canada.east", 
                  "SSE", "SSE.mpa" , "snowcrab")
bathymetry.db( p=p, DS="complete.redo" ) # finalise at diff resolutions 6 hrs ..
bathymetry.db( p=p, DS="baseline.redo" )  # coords of areas of interest ..filtering of areas and or depth to reduce file size, in planar coords only



### -----------------------------------------------------------------
# to recreate new polygons, run the following:
bathyclines.redo = FALSE
depths = c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 400, 500, 600, 700, 750, 800, 900,
             1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 )
if( bathyclines.redo ) {
  # note these polygons are created at the resolution specified in p$spatial.domain ..
  # which by default is very high ("canada.east.highres" = 0.5 km .. p$pres ).
  # For lower one specify an appropriate p$spatial.domain
  p = bio.bathymetry::bathymetry.parameters(resolution="canada.east") # reset to lower resolution 
  options(max.contour.segments=100000) # required if superhighres is being used
  plygn = isobath.db( p=p, DS="isobath.redo", depths=depths  )
}


### -----------------------------------------------------------------
p = bio.bathymetry::bathymetry.parameters(resolution="canada.east") # reset to lower resolution 
  
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

p = bio.bathymetry::bathymetry.parameters(resolution="canada.east.highres")
bathymetry.figures( p=p, vn="z", logyvar=TRUE )


p = bio.bathymetry::bathymetry.parameters(resolution="canada.east.superhighres")
bathymetry.figures( p=p, vn="z", logyvar=TRUE )


p = bio.bathymetry::bathymetry.parameters(resolution="snowcrab")
bathymetry.figures( p=p, vn="z", logyvar=TRUE )

