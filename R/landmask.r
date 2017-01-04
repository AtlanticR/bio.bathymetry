
  landmask = function( lonlat=NULL, db="worldHires", regions=c("Canada", "US"), proj4strvalue=NULL, return.value="not.land", tag="index", internal.projection=NULL, ... ) {
    #\\ Using the world coastline data base:
    #\\ return.value determines what is returned: "coast.lonat", "coast.polygon" (sp),
    #\\ or indices of"not.land", "land"

    if (is.null(internal.projection)) stop("internal.projection is required")

    require(maps)
    require(mapdata)
    require(maptools)
    require(rgdal)
    require(sp)

    fn = file.path( project.datadirectory("bio.bathymetry"), "landmask", paste( tag, db, paste0(regions, collapse=""), internal.projection, "rdata", sep="."))
    dir.create( dirname(fn), recursive=TRUE, showWarnings=FALSE  )

    if (is.null( lonlat)) {
      #\\ When lonlat is NULL, this is a flag to return a previous generated and saved version
      #\\ found in bio.data/bathymetry/landmask/ ...
      land = NULL
      if (file.exists( fn)) load(fn)
      if ( return.value=="test" )     return( land)
      if ( return.value=="not.land")  return( which ( is.na(land)) )
      if ( return.value=="land")      return( which ( !is.na(land)) )
    }

    if (is.null(proj4strvalue) ) proj4strvalue=CRS( sp::proj4string(lonlat) )
    if (is.null(proj4strvalue) ) proj4strvalue=CRS("+proj=longlat +datum=WGS84")

    coastline = maps::map( database=db, regions=regions, fill=TRUE, plot=FALSE, ...)
    if ( return.value=="coast.lonlat") return (coastline)

    coastlineSp = maptools::map2SpatialPolygons( coastline, IDs=coastline$names, proj4string=proj4strvalue  )
    if ( return.value=="coast.polygon") return (coastlineSp)

    land = over( SpatialPoints( lonlat, proj4strvalue ), coastlineSp )
    save( land, file=fn, compress=TRUE )
    return(fn)

  }


