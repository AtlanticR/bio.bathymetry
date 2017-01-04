
  bathymetry.db = function( p=NULL, DS=NULL, return.format="dataframe", varnames=NULL ) {

    datadir = project.datadirectory("bio.bathymetry", "data" )  # raw data
		dir.create( datadir, showWarnings=F, recursive=T )

    #\\ Note inverted convention: depths are positive valued
    #\\ i.e., negative valued for above sea level and positive valued for below sea level

    if ( DS=="gebco") {
      #library(RNetCDF)
      # request at: https://www.bodc.ac.uk/data/online_delivery/gebco/ [ jae.choi@dfo ] ... / gate.gate
      # extent: (WSEN) = -72,36,-45.,53
      # and saved as: bio.data/bathymetry/data/gebco.{xyz,nc}  # still waiting
      # and xz compressed
      fn = file.path( datadir, "bathymetry.gebco.rdata" )
      if (file.exists (fn) ) {
        load(fn)
        return(gebco)
      }

      fn_local = file.path( datadir, "gebco.xyz.xz") # xz compressed file
      nc = open.nc(bathy_fname)
      read.nc(nc)
      array(tmp$z, dim=tmp$dim)
      gebco = read.table( xzfile( fn_local ) )
      names(gebco) = c("lon", "lat", "z")
      gebco$z = - gebco$z
      # levelplot( log(z+ min(gebco$z) ))~lon+lat, gebco, aspect="iso")
      save( gebco, file=fn, compress=TRUE )
    }

    # --------------

    if ( DS=="etopo1") {
      # etopo1_bedrock.xyz ---> 1 min resolution
      # extent: (WSEN) = -72,36,-45.,53
      # download manually from:  http://maps.ngdc.noaa.gov/viewers/wcs-client/
      # and saved as: bio.data/bathymetry/data/etopo1_bedrock.xyz
      # and xz compressed
      fn = file.path( datadir, "bathymetry.etopo1.rdata" )
      if (file.exists (fn) ) {
        load(fn)
        return(etopo1)
      }
      fn_local = file.path( datadir, "etopo1_bedrock.xyz.xz") # xz compressed file
      etopo1 = read.table( xzfile( fn_local ) )
      names(etopo1) = c("lon", "lat", "z")
      etopo1$z = - etopo1$z
      # levelplot( log(z+ min(etopo1$z) ))~lon+lat, etopo1, aspect="iso")
      save( etopo1, file=fn, compress=TRUE )
    }

    # --------------

    if ( DS =="Greenlaw_DEM") {
      # DEM created 2014
      # GCS_WGS_1984, UTM_Zone_20N; spheroid:: 6378137.0, 298.257223563
      # 322624071 "grid points
      # 50 m  horizontal resolution
      # depth range: -5053.6 to 71.48 m
      fn = file.path( datadir, "bathymetry.greenlaw.rdata" )
      if (file.exists (fn) ) {
        load(fn)
        return(gdem)
      }

      require(rgdal)
      demfile.adf = file.path( datadir, "greenlaw_DEM", "mdem_50", "w001001.adf" )  # in ArcInfo adf format
      dem = new( "GDALReadOnlyDataset", demfile.adf )
      # gdem = asSGDF_GROD( dem, output.dim=dim(dem) ) # regrid to another dim
      # gdem = getRasterData(dem) # in matrix format
      gdem = getRasterTable(dem) # as a data frame
      names(gdem) = c("plon", "plat", "z")
      gdem = gdem[ is.finite( gdem$z ) , ]
      gdem$plon = gdem$plon / 1000
      gdem$plat = gdem$plat / 1000
      gdem = planar2lonlat( gdem, "utm20" )  # plon,plat in meters but crs for utm20 in km
      gdem = gdem[, c("lon", "lat", "z") ]
      save( gdem, file=project.datadirectory( "bio.bathymetry", "data", "bathymetry.greenlaw.rdata"), compress=TRUE )
    }

    # --------------

    if (  DS %in% c("z.lonlat.rawdata.redo", "z.lonlat.rawdata") ) {
			# raw data minimally modified all concatenated
      fn = file.path( datadir, "bathymetry.canada.east.lonlat.rawdata.rdata" )

      if (DS =="z.lonlat.rawdata" ) {
        load( fn )
        return( bathy )
      }

      print( "This is going to take a lot of RAM!")

			# this data was obtained from CHS via David Greenberg in 2004; range = -5467.020, 383.153; n=28,142,338
      fn_nwa = file.path( datadir, "nwa.chs15sec.xyz.xz") # xz compressed file
      chs15 = read.table( xzfile( fn_nwa ) )
      names(chs15) = c("lon", "lat", "z")
      # chs15 = chs15[ which( chs15$z < 1000 ) , ]
      chs15$z = - chs15$z

      # temporary break up of data to make it functional in smaller RAM systems
      chs1000_5000 = chs15[ which( chs15$z > 1000 ), ]
      u =  which(duplicated( chs1000_5000))
      if (length(u)>0) chs1000_5000 = chs1000_5000[-u,]

      chs0_1000 = chs15[ which( chs15$z <= 1000 ), ]
      u =  which(duplicated( chs0_1000 ))
      if (length(u)>0) chs0_1000 = chs0_1000[-u,]

      rm ( chs15); gc()

      fn0 = file.path( datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_0_1000.rdata" )
      fn1 = file.path( datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_1000_5000.rdata" )

      save ( chs0_1000, file=fn0 )
      save ( chs1000_5000, file=fn1 )

      rm ( chs0_1000, chs1000_5000 )
      gc()

      # pei = which( chs15$lon < -60.5 & chs15$lon > -64.5 & chs15$lat>45.5 & chs15$lat<48.5 )
      # levelplot( z~lon+lat, data=chs15[pei,] )


      # Michelle Greenlaw's DEM from 2014
      # range -3000 to 71.5 m; n=155,241,029 .. but mostly interpolated
      gdem = bathymetry.db( DS="Greenlaw_DEM" )
      gdem$z = - gdem$z

      # pei = which( gdem$lon < -60.5 & gdem$lon > -65 & gdem$lat>45.5 & gdem$lat<49 )
      # levelplot( z~I(round(lon,3))+I(round(lat,3)), data=gdem[pei,] )

      # bad boundaries in Greenlaw's gdem:
      # southern Gulf of St lawrence has edge effects
      bd1 = rbind( c( -62, 46.5 ),
                   c( -61, 47.5 ),
                   c( -65, 47.5 ),
                   c( -65, 46.2 ),
                   c( -64, 46.2 ),
                   c( -62, 46.5 ) )

      a = which( point.in.polygon( gdem$lon, gdem$lat, bd1[,1], bd1[,2] ) != 0 )
      gdem = gdem[- a,]

      # remove also the northern and eastern margins for edge effects
      gdem = gdem[ which( gdem$lat <  47.1) ,]
      gdem = gdem[ which( gdem$lon < -56.5) ,]

      fn0g = file.path( datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_0_1000_gdem.rdata" )
      fn1g = file.path( datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_1000_5000_gdem.rdata" )

      # temporary break up of data to make it functional in smaller RAM systems
      gdem1000_5000 = gdem[ which( gdem$z > 1000 ), ]
      # u =  which(duplicated( gdem1000_5000))
      # if (length(u)>0) gdem1000_5000 = gdem1000_5000[-u,]
      save ( gdem1000_5000, file=fn1g )
      rm( gdem1000_5000 ); gc()

      gdem0_1000 = gdem[ which( gdem$z <= 1000 ), ]
      # u =  which(duplicated( gdem0_1000 ))
      # if (length(u)>0) gdem0_1000 = gdem0_1000[-u,]
      save ( gdem0_1000, file=fn0g )
      rm( gdem0_1000 )
      rm( gdem) ;gc()
      gc()


      # chs and others above use chs depth convention: "-" is below sea level,
			# in snowcrab and groundfish convention "-" is above sea level
			# retain postive values at this stage to help contouring near coastlines

      bathy = bathymetry.db( DS="etopo1" )

      additional.data=c("snowcrab", "groundfish")

			if ( "snowcrab" %in% additional.data ) {
        # range from 23.8 to 408 m below sea level ... these have dropped the "-" for below sea level; n=5925 (in 2014)
        # bioLibrary( "bio.snowcrab")
        sc = bio.snowcrab::snowcrab.db("set.clean")[,c("lon", "lat", "z") ]
				sc = sc [ which (is.finite( rowSums( sc ) ) ) ,]
				j = which(duplicated(sc))
        if (length (j) > 0 ) sc = sc[-j,]
        bathy = rbind( bathy, sc )
			  # p = p0
        rm (sc); gc()

        #sc$lon = round(sc$lon,1)
        #sc$lat = round(sc$lat,1)
        # contourplot( z~lon+lat, sc, cuts=10, labels=F )
      }

      if ( "groundfish" %in% additional.data ) {
        # n=13031; range = 0 to 1054

        warning( "Should use bottom contact estimates as a priority ?" )
				gf = bio.groundfish::groundfish.db( "set.base" )[, c("lon","lat", "sdepth") ]
				gf = gf[ which( is.finite(rowSums(gf) ) ) ,]
        names(gf) = c("lon", "lat", "z")
				j = which(duplicated(gf))
        if (length (j) > 0 ) gf = gf[-j,]
 				bathy = rbind( bathy, gf )
        rm (gf); gc()

        #gf$lon = round(gf$lon,1)
        #gf$lat = round(gf$lat,1)
        #contourplot( z~lon+lat, gf, cuts=10, labels=F )

			}

      u =  which(duplicated( bathy ))
      if (length(u)>0) bathy = bathy[ -u, ]
      rm (u)

      bathy0 = bathy[ which(bathy$z <= 1000), ]
      bathy1 = bathy[ which(bathy$z  > 1000), ]
      rm(bathy)

      gc()

      load( fn0)
      bathy0 = rbind( bathy0, chs0_1000 )
      rm(chs0_1000) ;gc()

      load( fn0g )
      bathy0 = rbind( bathy0, gdem0_1000 )
      rm ( gdem0_1000 ) ;gc()

      u =  which(duplicated( bathy0 ))
      if (length(u)>0) bathy0 = bathy0[ -u, ]

      fn0b = file.path( datadir, "bathymetry.canada.east.lonlat.rawdata.temporary_0_1000_bathy.rdata" )
      save ( bathy0, file=fn0b )
      rm (bathy0); gc()

    # ---

      load( fn1 )
      bathy1 = rbind( bathy1, chs1000_5000 )
      rm (  chs1000_5000 ) ; gc()

      load( fn1g )
      bathy1 = rbind( bathy1, gdem1000_5000 )
      rm ( gdem1000_5000 ) ; gc()

      u =  which(duplicated( bathy1 ))
      if (length(u)>0) bathy1 = bathy1[ -u, ]
      rm (u)

      load( fn0b )

      bathy = rbind( bathy0, bathy1 )
      rm( bathy1, bathy0 ) ; gc()

      save( bathy, file=fn, compress=T )
      
      fn.bathymetry.xyz = file.path( project.datadirectory("bio.bathymetry"), "data", "bathymetry.canada.east.xyz" )
    
      if (exists("spatial.domain", p)) {
        # mechanism to override defaults, eg for other regions
        if ( p$spatial.domain %in% c("SSE", "snowcrab", "SSE.mpa", "canada.east", "canada.east.highres", "canada.east.superhighres" ) ) {
          fn.bathymetry.xyz = file.path( project.datadirectory("bio.bathymetry"), "data", "bathymetry.canada.east.xyz" )  # ascii
        }
      }

      fn.xz = xzfile( paste( fn.bathymetry.xyz, ".xz", sep="" ) )
      write.table( bathy, file=fn.xz, col.names=F, quote=F, row.names=F)
      system( paste( "xz",  fn.bathymetry.xyz ))  # compress for space

      file.remove( fn0)
      file.remove( fn1)
      file.remove( fn0g )
      file.remove( fn1g )
      file.remove( fn0b )

      return ( fn )
    }




		if ( DS %in% c( "Z.redo", "Z.lonlat", "Z.lonlat.grid", "Z.planar", "Z.planar.grid" ) ) {

			outdir = project.datadirectory("bio.bathymetry", "interpolated" )
			dir.create( outdir, showWarnings=F, recursive=T )

			fn.lonlat = file.path( outdir, paste( p$spatial.domain, "Z.interpolated.lonlat.rdata", sep=".") )
      fn.lonlat.grid = file.path( outdir, paste( p$spatial.domain, "Z.interpolated.lonlat.grid.rdata", sep=".") )
      fn.planar = file.path( outdir, paste(p$spatial.domain, "Z.interpolated.planar.rdata", sep=".") )
      fn.planar.grid = file.path( outdir, paste(p$spatial.domain, "Z.interpolated.planar.grid.rdata", sep=".") )

      if ( DS == "Z.lonlat" ) {
        load( fn.lonlat )
        return( Z )
      }
      if ( DS == "Z.lonlat.grid" ) {   # not used ... drop?
        load( fn.lonlat.grid )
        return( Z )
      }
      if ( DS == "Z.planar" ) {
        load( fn.planar )
        return( Z )
      }
      if ( DS == "Z.planar.grid" ) {    # used by map.substrate
        load( fn.planar.grid )
        return( Z )
      }

      Z0 = Z = bathymetry.db( p, DS="Z.gridded" )
      Z$lon = grid.internal( Z$lon, p$lons )
      Z$lat = grid.internal( Z$lat, p$lats )
      Z = block.spatial ( xyz=Z, function.block=block.mean )

      # create Z in lonlat xyz dataframe format
      save(Z, file=fn.lonlat, compress=T)

      # create Z in lonlat grid/matrix format
      Z = xyz2grid(Z, p$lons, p$lats)  # using R
      save(Z, file=fn.lonlat.grid, compress=T) # matrix/grid format

      # ---- convert to planar coords ...
      # must use the interpolated grid to get even higher resolution "data" with extrapolation on edges
      Z = lonlat2planar( Z0, proj.type= p$internal.projection )
      Z = Z[, c("plon", "plat", "z")]
      Z$plon = grid.internal( Z$plon, p$plons )
      Z$plat = grid.internal( Z$plat, p$plats )

      gc()
      Z = block.spatial ( xyz=Z, function.block=block.mean )

      # create Z in planar xyz format
      save( Z, file=fn.planar, compress=T )

      # create Z in matrix/grid format
      gc()
      Z = xyz2grid( Z, p$plons, p$plats)
      save( Z, file=fn.planar.grid, compress=T )

      return ("interpolated depths completed")

		}

		if ( DS %in% c( "dZ.redo", "dZ.lonlat", "dZ.lonlat.grid", "dZ.planar", "dZ.planar.grid" ) ) {

			outdir = file.path( project.datadirectory("bio.bathymetry"), "interpolated" )
			dir.create( outdir, showWarnings=F, recursive=T )

      fn.lonlat = file.path( outdir, paste( p$spatial.domain, "dZ.interpolated.lonlat.rdata", sep=".")  )
      fn.lonlat.grid = file.path( outdir, paste( p$spatial.domain, "dZ.interpolated.lonlat.grid.rdata", sep=".")  )
      fn.planar = file.path( outdir, paste( p$spatial.domain, "dZ.interpolated.planar.rdata", sep=".")  )
      fn.planar.grid = file.path( outdir, paste( p$spatial.domain, "dZ.interpolated.planar.grid.rdata", sep=".")  )

      if ( DS == "dZ.lonlat" ) {
        load( fn.lonlat )
        return( dZ )
      }
      if ( DS == "dZ.lonlat.grid" ) {
        load( fn.lonlat.grid )
        return( dZ )
      }
      if ( DS == "dZ.planar" ) {
        load( fn.planar )
        return( dZ )
      }
      if ( DS == "dZ.planar.grid" ) {
        load( fn.planar.grid )
        return( dZ )
      }

      dZ0 = bathymetry.db( p, DS="dZ.gridded" )
      dZ0$dZ = log( abs( dZ0$dZ ) )

      dZ = dZ0
      dZ$lon = grid.internal( dZ$lon, p$lons )
      dZ$lat = grid.internal( dZ$lat, p$lats )
      dZ = block.spatial ( xyz=dZ, function.block=block.mean )

      # create dZ in lonlat xyz dataframe format
      save(dZ, file=fn.lonlat, compress=T)

      # create dZ in lonlat grid/matrix format
      dZ = xyz2grid(dZ, p$lons, p$lats)
      save(dZ, file=fn.lonlat.grid, compress=T) # matrix/grid format


      # ---- convert to planar coords ...
      # must use the interpolated grid to get even higher resolution "data" with extrpolation on edges
      dZ = lonlat2planar( dZ0, proj.type= p$internal.projection )  # utm20, WGS84 (snowcrab geoid)
      dZ = dZ[, c("plon", "plat", "dZ")]
      dZ$plon = grid.internal( dZ$plon, p$plons )
      dZ$plat = grid.internal( dZ$plat, p$plats )

      gc()
      dZ = block.spatial ( xyz=dZ, function.block=block.mean )

      # create dZ in planar xyz format
      save( dZ, file=fn.planar, compress=T )

      # create dZ in matrix/grid format
      gc()
      dZ = xyz2grid( dZ, p$plons, p$plats)
      save( dZ, file=fn.planar.grid, compress=T )

      return ("interpolated files complete, load via another call for the saved files")
    }

    if ( DS %in% c( "ddZ.redo", "ddZ.lonlat", "ddZ.planar", "ddZ.lonlat.grid", "ddZ.planar.grid"  ) ) {
     	outdir = file.path( project.datadirectory("bio.bathymetry"), "interpolated" )
			dir.create( outdir, showWarnings=F, recursive=T )

      fn.lonlat = file.path( outdir,  paste( p$spatial.domain, "ddZ.interpolated.lonlat.rdata", sep=".")  )
      fn.lonlat.grid = file.path( outdir,  paste( p$spatial.domain, "ddZ.interpolated.lonlat.grid.rdata", sep=".")  )
      fn.planar = file.path( outdir,  paste( p$spatial.domain, "ddZ.interpolated.planar.rdata", sep=".")  )
      fn.planar.grid = file.path( outdir,  paste( p$spatial.domain, "ddZ.interpolated.planar.grid.rdata", sep=".") )

      if ( DS == "ddZ.lonlat" ) {
        load( fn.lonlat )
        return( ddZ )
      }
      if ( DS == "ddZ.lonlat.grid" ) {
        load( fn.lonlat.grid )
        return( ddZ )
      }
      if ( DS == "ddZ.planar" ) {
        load( fn.planar )
        return( ddZ )
      }
      if ( DS == "ddZ.planar.grid" ) {
        load( fn.planar.grid )
        return( ddZ )
      }

      ddZ0 = bathymetry.db( p, DS="ddZ.gridded" )
      ddZ0$ddZ = log( abs( ddZ0$ddZ ) )

      # ----- convert to lonlats blocked
      ddZ = ddZ0
      ddZ$lon = grid.internal( ddZ$lon, p$lons )
      ddZ$lat = grid.internal( ddZ$lat, p$lats )
      ddZ = block.spatial ( xyz=ddZ, function.block=block.mean )

      # lonlat xyz dataframe format
      save(ddZ, file=fn.lonlat, compress=T)

      ddZ = xyz2grid(ddZ, p$lons, p$lats)
      save(ddZ, file=fn.lonlat.grid, compress=T) # matrix/grid format

      # ---- convert to planar coords ...
      # must use the interpolated grid to get even higher resolution "data" with extrpolation on edges
      ddZ = ddZ0
      ddZ = lonlat2planar( ddZ0, proj.type= p$internal.projection )  # utm20, WGS84 (snowcrab geoid)
      ddZ = ddZ[, c("plon", "plat", "ddZ")]
      gc()
      ddZ$plon = grid.internal( ddZ$plon, p$plons )
      ddZ$plat = grid.internal( ddZ$plat, p$plats )

      ddZ = block.spatial ( xyz=ddZ, function.block=block.mean )

      # create ddZ in planar xyz format
      save( ddZ, file=fn.planar, compress=T )

      gc()
      ddZ = xyz2grid(ddZ, p$plons, p$plats)
      save( ddZ, file=fn.planar.grid, compress=T )

      return ("interpolated files complete, load via another call for the saved files")
    }

    # ------------

  
    if (DS %in% c("lookuptable.sse.snowcrab.redo", "lookuptable.sse.snowcrab" )) {
      #\\ DS="lookuptable.sse.snowcrab(.redo)" creates/returns a lookuptable for SSE -> snowcrab domains
      #\\   both share the same initial domains + resolutions and so it is faster to operate upon the indices
      fn = file.path( project.datadirectory("bio.bathymetry"), "interpolated", "sse.snowcrab.lookup.rdata")
      if (DS== "lookuptable.sse.snowcrab" ) {
        if (file.exists(fn)) load(fn)
        return(id)
      }
      zSSE = bathymetry.db ( p=spatial_parameters( type="SSE" ), DS="baseline" )
      zSSE$id.sse = 1:nrow(zSSE)

      zsc  = bathymetry.db ( p=spatial_parameters( type="snowcrab" ), DS="baseline" )
      zsc$id.sc = 1:nrow(zsc)

      z = merge( zSSE, zsc, by =c("plon", "plat"), all.x=T, all.y=T, sort=F )
      ii = which(is.finite(z$id.sc ) & is.finite(z$id.sse )  )
      if (length(ii) != nrow(zsc) ) stop("Error in sse-snowcrab lookup table size")
      id = sort( z$id.sse[ ii] )
      # oo= zSSE[id,]

      save( id, file=fn, compress=T )
      return(fn)
    }

    # ----------------

    if ( DS %in% c("lbm.inputs", "lbm.inputs.redo" )) {

      fn = file.path( datadir, paste( "bathymetry", "lbm.inputs", "rdata", sep=".") )
      if (DS =="lbm.inputs" ) {
        load( fn)
        return( hm )
      }
      
      print( "Warning: this needs a lot of RAM .. ~40GB depending upon resolution of discretization" )
      
      B = bathymetry.db ( p=p, DS="z.lonlat.rawdata" )  # 16 GB in RAM just to store!
      # B = B[ which(B$z > -100),]  # take part of the land to define coastline
      B = lonlat2planar( B, proj.type=p$internal.projection )
      B$lon = NULL
      B$lat = NULL
      gc()
      
      # if using kernel density/FFT.. might as well take the gird here as it requires  gridded data 
      B$plon = grid.internal( B$plon, p$plons )
      B$plat = grid.internal( B$plat, p$plats )
      B = block.spatial ( xyz=B[,c("plon", "plat", "z")], function.block=block.mean )

#       B = B[, c("plon", "plat", "z")]

      BP = list( LOCS = expand.grid( p$plons, p$plats, KEEP.OUT.ATTRS=FALSE) )
      names( BP$LOCS ) = c("plon", "plat")

      hm = list( input=B, output=BP )
      save( hm, file=fn, compress=TRUE)
    }

    # ----------------

    if ( DS == "landmasks.create" ) {
      # on resolution of predictions
      pps  =  expand.grid( plons=p$plons, plats=p$plats, KEEP.OUT.ATTRS=FALSE)
      V = SpatialPoints( planar2lonlat( pps, proj.type=p$internal.crs )[, c("lon", "lat" )], CRS("+proj=longlat +datum=WGS84") )
      landmask( lonlat=V, db="worldHires", regions=c("Canada", "US"), ylim=c(36,53), xlim=c(-72,-45), tag="predictions" ,internal.projection=p$internal.projection)

      # on resolution of statistics
      V = as.data.frame(p$ff$Sloc[])
      names(V) = c("plon", "plat")
      # V = data.frame( cbind(plon=Sloc[,1], plat=Sloc[,2]) )
      V = SpatialPoints( planar2lonlat( V, proj.type=p$internal.crs )[, c("lon", "lat" )], CRS("+proj=longlat +datum=WGS84") )
      landmask( lonlat=V, db="worldHires",regions=c("Canada", "US"), ylim=c(36,53), xlim=c(-72,-45), tag="statistics" ,internal.projection=p$internal.projection)
    }

    #-------------------------

    if ( DS %in% c("lbm.finalize.redo", "lbm.finalize" )) {
      #// bathymetry( p, DS="lbm.finalize(.redo)" return/create the
      #//   lbm interpolated method formatted and finalised for production use
      fn = file.path(  project.datadirectory("bio.bathymetry"), "interpolated",
        paste( "bathymetry", "lbm", "finalized", p$spatial.domain, "rdata", sep=".") )
      if (DS =="lbm.finalize" ) {
        B = NULL
        if ( file.exists ( fn) ) load( fn)
        return( B )
      }

      nr = p$nplons
      nc = p$nplats

      # data prediction grid
      B = expand.grid( plon=p$plons, plat=p$plats, KEEP.OUT.ATTRS=FALSE) 

      Bmean = lbm_db( p=p, DS="lbm.prediction", ret="mean" )
      Bsd = lbm_db( p=p, DS="lbm.prediction", ret="sd" )
      B = cbind(B, Bmean, Bsd)
      rm (Bmean, Bsd); gc()
      names(B) = c( "plon", "plat", "z", "z.sd") # really Z.mean but for historical compatibility "z"

      # # remove land
      # oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="land", tag="predictions" )
      # B$z[oc] = NA
      # B$z.sd[oc] = NA

      Bmn = matrix( data=B$z, nrow=nr, ncol=nc )  # means

      # first order central differences but the central term drops out:
      # diffr = ( ( Bmn[ 1:(nr-2), ] - Bmn[ 2:(nr-1), ] ) + ( Bmn[ 2:(nr-1), ] - Bmn[ 3:nr, ] ) ) / 2
      # diffc = ( ( Bmn[ ,1:(nc-2) ] - Bmn[ ,2:(nc-1) ] ) + ( Bmn[ ,2:(nc-1) ] - Bmn[ ,3:nc ] ) ) / 2
      diffr =  Bmn[ 1:(nr-2), ] - Bmn[ 3:nr, ]
      diffc =  Bmn[ ,1:(nc-2) ] - Bmn[ ,3:nc ]
      rm (Bmn); gc()

      dZ = ( diffr[ ,2:(nc-1) ] + diffc[ 2:(nr-1), ] ) / 2
      dZ = rbind( dZ[1,], dZ, dZ[nrow(dZ)] )  # top and last rows are copies .. dummy value to keep dim correct
      dZ = cbind( dZ[,1], dZ, dZ[,ncol(dZ)] )

      B$dZ =  abs(c(dZ))

      # gradients
      ddiffr =  dZ[ 1:(nr-2), ] - dZ[ 3:nr, ]
      ddiffc =  dZ[ ,1:(nc-2) ] - dZ[ ,3:nc ]
      rm( dZ ); gc()

      ddZ = ( ddiffr[ ,2:(nc-1) ] + ddiffc[ 2:(nr-1), ] ) / 2
      ddZ = rbind( ddZ[1,], ddZ, ddZ[nrow(ddZ)] )  # top and last rows are copies .. dummy value to keep dim correct
      ddZ = cbind( ddZ[,1], ddZ, ddZ[,ncol(ddZ)] )
      B$ddZ = abs(c(ddZ))

      # merge into statistics
      BS = lbm_db( p=p, DS="stats.to.prediction.grid" )
      B = cbind( B, BS )
      rm (BS); gc()

      # names(B) = c( names(B), p$statvars )

      save( B, file=fn, compress=TRUE)
      return(fn)

      if (0) {
        aoi = which( B$z > 5 & B$z < 3000 )
        levelplot( log(z) ~ plon + plat, B[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( log(phi) ~ plon + plat, B[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( log(range) ~ plon + plat, B[aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )

      }
    }

    # ------------


    if ( DS %in% c( "complete", "complete.redo") ) {
      #// everything, including land 

      Z = NULL

      if ( DS %in% c( "complete") ) {

        domain = NULL
        if ( is.null(domain)) {
          if ( exists("spatial.domain", p)) {
            domain = p$spatial.domain
          } else if ( exists( "new.grids", p) ) { # over-rides p$spatial domain
            if( length( p$new.grids )== 1 ) {
              domain = p$new.grids
        } } }

        fn = file.path( project.datadirectory("bio.bathymetry", "interpolated"),
          paste( "bathymetry", "complete", domain, "rdata", sep=".") )
        if ( file.exists ( fn) ) load( fn)
        return( Z )

      }

      Z = list()

      p0 = p  # the originating parameters
      Z0 = bathymetry.db( p=p0, DS="lbm.finalize" )
      coordinates( Z0 ) = ~ plon + plat
      crs(Z0) = crs( p0$interal.crs )
      
      grids = unique( c( p$spatial.domain, p$new.grids ))

      for (gr in grids ) {
        print(gr)
        p1 = spatial_parameters( type=gr )
        for (vn in names(Z0)) {
          Z[[vn]] = raster::projectRaster(
            from = raster::rasterize( Z0, bio.spacetime::spatial_parameters_to_raster(p0), field=vn, fun=mean),
            to   = bio.spacetime::spatial_parameters_to_raster( p1) )
        }
        Z = as( brick(Z), "SpatialPointsDataFrame" )
        Z = as.data.frame(Z)
        u = names(Z)
        names(Z)[ which( u=="x") ] = "plon"
        names(Z)[ which( u=="y") ] = "plat"
        fn = file.path( project.datadirectory("bio.bathymetry", "interpolated"),
          paste( "bathymetry", "complete", p1$spatial.domain, "rdata", sep=".") )
        save (Z, file=fn, compress=TRUE)
      }


if (0) {
  # -- incomplete TODO .. need a layer to covert everything to the same projection first before warping 
  # and then return it to desired output projection  
  
      p0 = p  # the originating parameters

      p0$wgts = fields::setup.image.smooth( 
        nrow=p0$nplons, ncol=p0$nplats, dx=p0$pres, dy=p0$pres,
        theta=p0$pres, xwidth=4*p0$pres, ywidth=4*p0$pres )
      
      Z0 = bathymetry.db( p=p0, DS="lbm.finalize" )
      coords = c("plon","plat")
      varnames = setdiff( names(Z0), coords )  
      #using fields

      grids = unique( c( p$spatial.domain, p$new.grids ))
      Z0_i = as.matrix( round( cbind( 
        ( Z0$plon-p0$plons[1])/p0$pres + 1, (Z0$plat-p0$plats[1])/p0$pres + 1
      ) ) ) 

      for (gr in grids ) {
        print(gr)
        p1 = spatial_parameters( type=gr ) #target projection
        Z = expand.grid( plon=p1$plons, plat=p1$plats, KEEP.OUT.ATTRS=FALSE )

        for (vn in varnames) {
          M = fields::as.image( Z0[[vn]], ind=Z0_i, na.rm=TRUE, nx=p0$nplons, ny=p0$nplats )
          Z1 = fields::interp.surface( M, loc=Z[,coords] )
          Z[[vn]] = c(Z1)
          ii = which( !is.finite( Z[[vn]] ) )
          if ( length( ii) > 0 ) {
            # try again ..
            Z2 = fields::interp.surface( M, loc=Z[ii, coords] )
            Z[[vn]][ii] = Z2[ii]  
          }
          ii = which( !is.finite( Z[[vn]] ) )
          if ( length( ii) > 0 ) {
            Z3 =  fields::image.smooth( M, dx=p0$pres, dy=p0$pres, wght=p0$wght )$z
            Z[[vn]][ii] = Z3[ii]
          }
        }

        fn = file.path( project.datadirectory("bio.bathymetry", "interpolated"),
          paste( "bathymetry", "complete", p1$spatial.domain, "rdata", sep=".") )
        save (Z, file=fn, compress=TRUE)
      }


}



      return ( "Completed subsets" )
    }



    # ------------

    if (DS %in% c("baseline", "baseline.redo") ) {
      # form prediction surface in planar coords over the ocean

      if ( DS=="baseline" ) {
        outfile =  file.path( project.datadirectory("bio.bathymetry"), "interpolated",
          paste( p$spatial.domain, "baseline.interpolated.rdata" , sep=".") )
        Z = NULL
        load( outfile )
        if (is.null(varnames)) varnames =c("plon", "plat")
        Znames = names(Z)
        varnames = intersect( Znames, varnames )
        Z = Z[ , varnames]
        return (Z)
      }

      for (domain in p$new.grids) {
        pn = spatial_parameters( type=domain )
        if ( pn$spatial.domain == "snowcrab" ) {
          # NOTE::: snowcrab baseline == SSE baseline, except it is a subset so begin with the SSE conditions
          pn = spatial_parameters( type="SSE", p=pn ) 
        }
        Z = bathymetry.db( p=pn, DS="complete"  )
        Z = filter.bathymetry( DS=domain, Z=Z )

        outfile =  file.path( project.datadirectory("bio.bathymetry"), "interpolated",
          paste( domain, "baseline.interpolated.rdata" , sep=".") )

        save (Z, file=outfile, compress=T )
        print( outfile )
      }
      # require (lattice); levelplot( z~plon+plat, data=Z, aspect="iso")
      return( "completed" )
    }



  }  # end bathymetry.db








