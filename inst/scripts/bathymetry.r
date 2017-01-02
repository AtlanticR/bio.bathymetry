
# Bathymetry 
# warning: this will take weeks as it is an iterative process


if ( basedata.redo ) {
  p = bio.bathymetry::bathymetry.parameters() 
  bathymetry.db( p=p, DS="z.lonlat.rawdata.redo" ) # needs about 42 GB RAM, JC 2015
  bathymetry.db( p=p, DS="landmasks.create" ) # re-run only if default resolution .. v. slow ... currently using sp::over ... replace with point.in.polygon (TODO)
  bathymetry.db( p=p, DS="bathymetry.lbm.redo" )  # Warning: req ~ 15 min, 40 GB RAM (2015, Jae) data to model (with covariates if any)
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


p = lbm( p=p, DATA='bathymetry.db( p=p, DS="bathymetry.lbm" )' )   
if (restarting) {
  lbm_db(p=p, DS="statistics.status.reset" )
  p = lbm( p=p, continue=TRUE ) 
}


# bring together stats and predictions and any other required computations: slope and curvature
bathymetry.db( p=p, DS="bathymetry.lbm.finalize.redo" ) 
# B = bathymetry( p=p, DS="bathymetry.lbm.finalize" )     # to see the assimilated data:


# as the interpolation process is so expensive, regrid/upscale/downscale based off the above run
# if you want more, will need to add to the list and modify the selection criteria
p$new.grids = c( "canada.east.ultrahighres", "canada.east.highres", "canada.east", "SSE", "SSE.mpa" , "snowcrab")
bathymetry.db( p=p, DS="complete.redo", grids.new=p$new.grids )
bathymetry.db( p=p, DS="baseline.redo" )   # filtering of areas and or depth to reduce file size, in planar coords only


# "snowcrab" subsets do exist but are simple subsets of SSE
# so only the lookuptable below is all that is important as far as bathymetry is concerned
# both share the same initial domains + resolutions
p = bio.bathymetry::bathymetry.parameters(resolution="snowcrab") # reset to defaults
bathymetry.db( p=p, DS="lookuptable.sse.snowcrab.redo" ) # indices to map SSE to snowcrab


### -----------------------------------------------------------------
# to recreate new polygons, run the following:
bathyclines.redo = FALSE
depths = c( 0, 10, 20, 50, 75, 100, 200, 250, 300, 400, 500, 600, 700, 750, 800, 900,
             1000, 1200, 1250, 1400, 1500, 1750, 2000, 2500, 3000, 4000, 5000 )
if( bathyclines.redo ) {
  # note these polygons are created at the resolution specified in p$spatial.domain ..
  # which by default is very high ("canada.east.highres" = 0.5 km .. p$pres ).
  # For lower one specify an appropriate p$spatial.domain
  p = bio.bathymetry::bathymetry.parameters() # reset to defaults
  plygn = isobath.db( p=p, DS="isobath.redo", depths=depths  )
}


### -----------------------------------------------------------------
p = bio.bathymetry::bathymetry.parameters() # reset to defaults
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
p = bio.bathymetry::bathymetry.parameters(resolution="canada.east.highres") # reset to defaults
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



