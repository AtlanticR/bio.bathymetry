
bathymetry.parameters = function(DS, p=NULL, resolution="canada.east.highres", nc=1 ) {

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
    p$bathymetry.bigmemory.reset = FALSE
    # cluster definition
    p$clusters = rep( "localhost", nc )

    return(p)
  }

}

