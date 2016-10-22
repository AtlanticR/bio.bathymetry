
  bathymetry.db = function( p=NULL, DS=NULL, additional.data=c("snowcrab", "groundfish"), return.format="dataframe" ) {

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
      nc <- open.nc(bathy_fname)
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

      fn.xz = xzfile( paste( p$bathymetry.xyz, ".xz", sep="" ) )
      write.table( bathy, file=fn.xz, col.names=F, quote=F, row.names=F)
      system( paste( "xz",  p$bathymetry.xyz ))  # compress for space

      file.remove( fn0)
      file.remove( fn1)
      file.remove( fn0g )
      file.remove( fn1g )
      file.remove( fn0b )

      return ( fn )
    }


    # --------------

    if ( DS == "gmt.retired" ) {
      # *** THIS is deprecated. ***
      # *** THIS is deprecated. ***
      # *** THIS is deprecated. ***

      # It is stored for historical purposes or as a fast method if you do not have access to
      # large computational resources

      p=list()
      p$libs = bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry" )
      p$libs = c( p$libs, RLibrary( "rgdal", "lattice", "parallel" ) )

      # ------------------
      # glue all data sources (spherical coords)
      # ... right now this is about 17 GB in size when expanded .... SLOW ....
      # and it takes about 52+ GB RAM (due to addition of Greenlaw's DEM )
      # run on servers only unless your machine can handle it
      redo.bathymetry.rawdata = FALSE
      if ( redo.bathymetry.rawdata ) {
        p = spacetime_parameters( type="canada.east", p=p )
        bathymetry.db ( p, DS="z.lonlat.rawdata.redo", additional.data=c("snowcrab", "groundfish") )
     }

      # ------------------
      # GMT-based methods:
      # NOTE: too many clusters will overload the system as data files are large ~(11GB RAM required to block)
      # for the high resolution maps .. the temporary files can be created deleted/overwritten files
      # in the temporary drives
      redo.isobaths = FALSE
      if (redo.isobaths) {
        area = c( "snowcrab", "SSE", "ecnasap", "canada.east" )
        for (sp in area) {
          p$spatial.domain = sp
          p = spacetime_parameters( p=p )
          p = gmt.parameters(p)  # interpolation parameters ... currently using GMT to interpolate bathymetry
          # override defaults in gmt.parameters as additional ones are used by other systems including lattice
          p$isobaths = c( 0, seq(50, 450, by=50), seq( 500, 1000, by=100 )  ) #override defaults
          p = make.list( list( depths=p$isobaths ), Y=p )
          p$clusters = rep( "localhost", 1 )
          #isobath.db( p=p, DS="gmt.redo" )
          parallel.run( isobath.db,  p=p, DS="gmt.redo" )
        }
      }

      # ------------------
      # intermediary base maps with location definitions, annotations and isobaths ... to speed up PS map production .. only for GMT maps
      redo.basemap.gmt = FALSE
      if ( redo.basemap.gmt ) {
        area = c( "snowcrab", "SSE", "ecnasap", "canada.east" )
        for (sp in area) {
          p$spatial.domain = sp
          p = spacetime_parameters( p=p )
          p = gmt.parameters(p)  # interpolation parameters ... currently using GMT's isobaths whcih are specified in gmt.parameters
          # or if you want to override the isobaths plotted define them here (but make sure they were created in the previous step)
          # p$isobaths = c( seq(50, 450, by=100)  )
          gmt.basemap(p)
        }
      }

      # ------------------
      # prepare finalised bathymetry data for use in ecomod
      complete.bathymetry.db = FALSE
      areas = c( "canada.east", "SSE" ) # only two are currently used
      for ( sp in areas ) {
        p = spacetime_parameters( type=sp, p=p )
        p = gmt.parameters(p)
        bathymetry.db ( p, DS="prepare.intermediate.files.for.dZ.ddZ" )  # uses GMT's math functions ...
        bathymetry.db ( p, DS="Z.redo" )
        bathymetry.db ( p, DS="dZ.redo" )
        bathymetry.db ( p, DS="ddZ.redo" )
        bathymetry.db ( p, DS="baseline.redo" ) # additional filtering of areas and or depth to reduce file size
        bathymetry.db ( p, DS="complete.redo" ) # glue all together
      }

      # ------------------
      # "snowcrab" subsets do exist but are simple subsets of SSE
      # so only the lookuptable below is all that is important as far as bathymetry is concerned
      # both share the same initial domains + resolutions
      p = spacetime_parameters( type="snowcrab", p=p )
      bathymetry.db ( p, DS="baseline.redo" ) # additional filtering of areas and or depth to reduce file size
      bathymetry.db( DS="lookuptable.sse.snowcrab.redo" )

      # ------------------
      ## a few lattice-based maps: for SSE only right now
      p = spacetime_parameters( type="SSE" )
      x = bathymetry.db ( p, DS="baseline" )

      snowcrab.area=F
      if (snowcrab.area) {
        # this is used below
        sc = intersect(
            which( x$plon< 990 & x$plon > 220  & x$plat< 5270 & x$plat > 4675 ) ,
            filter.region.polygon( x[, c("plon", "plat") ], "cfaall", planar=T)
        )
        x = x[sc,]
      }

      x$z =log( x$z )

      outdir = file.path(project.datadirectory("bio.bathymetry","maps"), p$spatial.domain)

      dr = quantile( x$z, probs=c(0.005, 0.995))
      datarange = seq(dr[1], dr[2], length.out=100)
      cols = color.code( "blue.black", datarange )
      outfn = "depth"
      annot = "ln ( Depth; m )"
      map( xyz=x[,c("plon", "plat", "z")], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
        fn=outfn, loc=outdir, at=datarange , col.regions=cols, spatial.domain=p$spatial.domain )


      x = bathymetry.db ( p, DS="dZ.planar" )
      if (snowcrab.area) x = x[sc,]
      dr = quantile( x$dZ, probs=c(0.005, 0.995))
      datarange = seq(dr[1], dr[2], length.out=100)
      cols = color.code( "blue.black", datarange )
      outfn = "slope"
      annot = "ln ( Slope; m/m )"
      map( xyz=x[ ,c("plon", "plat", "dZ")], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
        fn=outfn, loc=outdir, at=datarange , col.regions=cols , spatial.domain=p$spatial.domain )


      x = bathymetry.db ( p, DS="ddZ.planar" )
      if (snowcrab.area) x = x[sc,]
      dr = quantile( x$ddZ, probs=c(0.005, 0.995))
      datarange = seq(dr[1], dr[2], length.out=100)
      cols = color.code( "blue.black", datarange )
      outfn = "curvature"
      annot = "ln ( Curvature; m/m/m )"
      map( xyz=x[,c("plon", "plat", "ddZ")], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
        fn=outfn, loc=outdir, at=datarange , col.regions=cols, spatial.domain=p$spatial.domain )

      # Googleearth overlays:

      # common functions:
      bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry")

      cmd = function(x, ...) { system(paste(x, ...)) }

      tmpdir = tempdir()

      # 1.  map of bathymetry contours : colour background and lines

      gmt = list()
      gmt$out = file.path(project.datadirectory("bio.bathymetry"), "ss.bathymetry.colour.platecarre.ps" )
      gmt$outputs = c( "colourscale", "colourcontour", "bathymetry.redo" )
      gmt$region=" -R-72/-52/40/50"
      gmt$gmtproj="-JQ-62/6.5i"  # Cylindrical equidistant (Plate Carre)  the default required by Google Earth
      gmt$resolution="-I0.25m" #  resolution
      gmt$inp = file.path(project.datadirectory("bio.bathymetry"), "data", "bathymetry.canada.east.xyz" )
      gmt$tension ="-T0.35"
      gmt$cpt = "-Cjet -T-350/350/10 -Z -D"
      gmt$incscale = "-B50"
      gmt$font = "14p"
      gmt$annot.base = "-63.7 47.25 24 0 Helvetica LT"  # lon0, lat0, fontsize, angle, font, justification
      gmt$annot.text = "Bathymetry (m)"
      gmt$scale.location = "-D4.75i/0.75i/2.6i/0.16ih" # alternate: "-D4.5i/0.8i/2.5i/.25ih"
      # gmt$scale.location = "-D4.75i/0.75i/1i/0.1ih" # alternate: "-D4.5i/0.8i/2.5i/.25ih"

      setwd( dirname( gmt$out ) )
      gmt.map.simple(gmt)
      cmd( "ps2raster", gmt$out, "-V -Au -P -E300 -Tg -S" )

      #  cmd( "ps2raster", gmt$out, "-V -Au -P -E300 -TG -Qg2 -S" )

     # 2. map point-kriged data .. used for overlaying onto google-earth to id new stations ---------
     # fishable biomass density

      gmt = bio.snowcrab::load.environment()

      gmt$out = file.path( project.datadirectory("bio.snowcrab"), "maps", "googleearth", "R0.platecarre.ps" )
      gmt$outputs = c( "colourscale", "colourcontour" )
      gmt$region=" -R-72/-52/40/50"
      gmt$gmtproj="-JQ-62/6.5i"  # Cylindrical equidistant (Plate Carre)  the default required by Google Earth
      gmt$resolution="-I0.25m" #  resolution
      gmt$inp = NULL
      gmt$tension ="-T1"

      gmt$cpt = "-Cseis -T-3/2/0.1 -Z -I -D"
      gmt$incscale = "-B2"
      gmt$font = "14p"
      gmt$annot.base = "-63.7 47.25 24 0 Helvetica LT"  # lon0, lat0, fontsize, angle, font, justification
      gmt$annot.text = "Fishable biomass (log10; kg/km2)"
      gmt$scale.location = "-D4.25i/0.75i/0.75i/0.075ih" # alternate: "-D4.5i/0.8i/2.5i/.25ih"

      gmt$bathy.xyz = file.path( project.datadirectory("bio.bathymetry"), "data", "bathymetry.canada.east.xyz")
      gmt$bathy.grid = file.path(tmpdir, make.random.string(".gmt.grid"))

      PS = kriging.db( DS="UK.point.PS", p=list(v="R0.mass", y=2008, r="cfaall", transgaussian.kriging=T)  )
      PS$R0.mass.pred = log10( PS$R0.mass.pred )
      PS = planar2lonlat( PS, proj.type="utm20" )

      gmt$dat = PS[, c("lon", "lat", "R0.mass.pred" )]

      # determine colour ranges:
      hist( gmt$dat[,3]) # log scale ssems to work best: -3 to 2

      gmt.map.simple(gmt)

      setwd( dirname( gmt$out ) ) # ps2raster has difficulties with different directories

      cmd( "ps2raster", basename(gmt$out), "-V -Au -P -E300 -Tg -S" )
      #  cmd( "ps2raster", gmt$out, "-V -Au -P -E300 -TG -Qg2 -S" )

      # then use Gimp or something like it to replace white with transparent

     # 3. map point-kriged data .. used for overlaying onto google-earth to id new stations ---------
     # SD in fishable biomass

      gmt = bio.snowcrab::load.environment()
      gmt$out = file.path( project.datadirectory("bio.snowcrab"), "maps", "googleearth", "R0.sd.platecarre.ps" )
      gmt$outputs = c( "colourscale", "colourcontour" )
      gmt$region=" -R-72/-52/40/50"
      gmt$gmtproj="-JQ-62/6.5i"  # Cylindrical equidistant (Plate Carre)  the default required by Google Earth
      gmt$resolution="-I0.25m" #  resolution
      gmt$inp = NULL
      gmt$tension ="-T1"
      gmt$bathy.xyz = file.path( project.datadirectory("bio.bathymetry"), "data", "bathymetry.canada.east.xyz")
      gmt$bathy.grid = file.path(tmpdir, make.random.string(".gmt.grid"))

      PS = kriging.db( DS="UK.point.PS", p=list(v="R0.mass", y=2008, r="cfaall", transgaussian.kriging=T)  )
      PS = planar2lonlat( PS, proj.type="utm20" )
      PS$R0.mass.sd = sqrt( PS$R0.mass.var )
      gmt$dat = PS[, c("lon", "lat", "R0.mass.sd" )]
      gmt$dat = gmt$dat[ which( gmt$dat[,3] <=14 ) ,]
      # determine colour ranges:
      hist( gmt$dat[,3] )
      hist( log10( gmt$dat[,3] ) ) # log scale ssems to work best: -3 to 2

      gmt$cpt = "-Cseis -T1/12/0.1 -Z -I -D"
      gmt$incscale = "-B2"
      gmt$font = "14p"
      gmt$annot.base = "-63.7 47.25 24 0 Helvetica LT"  # lon0, lat0, fontsize, angle, font, justification
      gmt$annot.text = "Fishable biomass SD (log10; kg/km2)"
      gmt$scale.location = "-D4.25i/0.75i/0.75i/0.075ih" # alternate: "-D4.5i/0.8i/2.5i/.25ih"

      gmt.map.simple(gmt)

      setwd( dirname( gmt$out ) ) # ps2raster has difficulties with different directories
      #  cmd( "ps2raster", gmt$out, "-V -Au -P -E300 -TG -Qg2 -S" )
      cmd( "ps2raster", basename(gmt$out), "-V -Au -P -E300 -Tg -S" )

      ### END GMT-based methods
    }

    # --------------

    if ( DS %in% c("prepare.intermediate.files.for.dZ.ddZ", "Z.gridded", "dZ.gridded", "ddZ.gridded" ) ) {

			tmpdir  = tempdir()
			outdir = project.datadirectory("bio.bathymetry", "interpolated" )
			dir.create( outdir, showWarnings=F, recursive=T )

			fn.interp.z = file.path( outdir,  paste(  p$spatial.domain, "Z.interpolated.xyz", sep=".") )
			fn.interp.dz = file.path( outdir,  paste( p$spatial.domain, "dZ.interpolated.xyz", sep=".") )
			fn.interp.ddz = file.path( outdir,  paste( p$spatial.domain, "ddZ.interpolated.xyz", sep=".") )

			if (DS=="Z.gridded") {
				Z = read.table( fn.interp.z )
				names( Z ) = c("lon", "lat", "z")
				return( Z )
			}
			if (DS=="dZ.gridded") {
				dZ = read.table( fn.interp.dz )
				names( dZ ) = c("lon", "lat", "dZ")
				return( dZ )
			}
			if (DS=="ddZ.gridded") {
				ddZ = read.table( fn.interp.ddz )
				names( ddZ ) = c("lon", "lat", "ddZ")
				return( ddZ )
			}

			append = "-O -K"
			b.res = "-I10s"  # use full resolution of bathymetry data
			bathy.tension = "-T0.75"  # large steepness :: 0.35+ for steep; 0.25 for smooth
			blocked = file.path(tmpdir, make.random.string(".gmt.blocked"))
			grids  = file.path(tmpdir, make.random.string( ".gmt.depths"))
			z.dds = file.path(tmpdir, make.random.string( ".gmt.z.dds"))
			z.d2ds2 = file.path(tmpdir, make.random.string( ".gmt.z.d2ds2"))

      if ( !file.exists( p$bathymetry.bin )) {
        # a GMT binary file of bathymetry .. currently, only the "canada.east" domain
        # is all that is required/available
        fnxz = paste(p$bathymetry.xyz, ".xz", sep="")
        cmd( "xz --decompress", fnxz )
        cmd( "gmtconvert -bo", p$bathymetry.xyz, ">", p$bathymetry.bin )
        cmd( "xz --compress", p$bathymetry.xyz )
      }

			cmd( "blockmean", p$bathymetry.bin, "-bi3 -bo", p$region, b.res, ">", blocked )
			cmd( "surface", blocked, "-bi3", p$region, b.res, bathy.tension, paste("-G", grids, sep="") )
			cmd( "grdmath -M", grids, "DDX ABS", grids, "DDY ABS ADD 0.5 MUL =", z.dds )
			cmd( "grdmath -M -N", grids, "CURV =", z.d2ds2 )
			cmd( "grd2xyz", grids, ">", fn.interp.z )  # the scalar in meter
			cmd( "grd2xyz", z.dds, ">", fn.interp.dz )  # the scalar in meter / meter
			cmd( "grd2xyz", z.d2ds2, ">", fn.interp.ddz )  # the scalar m / m^2

			remove.files( c(blocked, grids, z.dds, z.d2ds2  ) )
			return ("intermediate files completed")
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

    if (DS %in% c("baseline.gmt", "baseline.gmt.redo") ) {
      #\\ form prediction surface in planar coords for SS snowcrab area
      #\\ deprecated .. use "baseline" (below)
      outfile =  file.path( project.datadirectory("bio.bathymetry"), "interpolated",
          paste( p$spatial.domain, "baseline.interpolated.gmt.rdata" , sep=".") )

      if ( DS=="baseline.gmt" ) {
        load( outfile )
        return (Z)
      }

      if ( p$spatial.domain == "snowcrab" ) {
        # NOTE::: snowcrab baseline == SSE baseline, except it is a subset so begin with the SSE conditions
        Z = bathymetry.db( p=spacetime_parameters( type="SSE", p=p ), DS="baseline.gmt" )
      } else {
        Z = bathymetry.db( p=p, DS="Z.planar" )
      }

      Z = filter.bathymetry( DS=p$spatial.domain, Z=Z )
      save (Z, file=outfile, compress=T )
			return( paste( "Baseline data file completed:", outfile )  )
      # require (lattice); levelplot( z~plon+plat, data=Z, aspect="iso")
    }

    # ------------

    if (DS %in% c("baseline", "baseline.redo") ) {
      # form prediction surface in planar coords
      outfile =  file.path( project.datadirectory("bio.bathymetry"), "interpolated",
          paste( p$spatial.domain, "baseline.interpolated.rdata" , sep=".") )

      if ( DS=="baseline" ) {
        load( outfile )
        return (Z)
      }

      for (domain in p$new.grids) {
        pn = spacetime_parameters( type=domain )
        if ( pn$spatial.domain == "snowcrab" ) {
          # NOTE::: snowcrab baseline == SSE baseline, except it is a subset so begin with the SSE conditions
          Z = bathymetry.db( p=spacetime_parameters( type="SSE", p=pn ), DS="complete", return.format = "dataframe.filtered"  )
        } else {
          Z = bathymetry.db( p=pn , DS="complete", return.format = "dataframe.filtered"  )
        }
        names0 = names( Z)
        Z = as.data.frame(Z)
        names(Z) = c( names0, "plon", "plat")
        Z = Z[, c("plon", "plat", "z")]
        save (Z, file=outfile, compress=T )
        print( outfile )
      }
      # require (lattice); levelplot( z~plon+plat, data=Z, aspect="iso")
			return( "completed" )
    }


    # --------------

    if (DS %in% c( "complete.gmt", "complete.gmt.redo") ) {
      #\\ DS="complete(.redo)" creates or returns the prediction surface in planar coords for SS snowcrab area
      #\\  derived from GMT-bsaed methods, deprecated ... use "complete" (below)

      # following methods are now deprecated and here on for backaward compatibility
      # to access use DS="complete.gmt"

      outfile =  file.path( project.datadirectory("bio.bathymetry"), "interpolated", paste( p$spatial.domain, "complete.rdata" , sep=".") )
      if (p$spatial.domain == "snowcrab" ) outfile=gsub( p$spatial.domain, "SSE", outfile )

      if ( DS=="complete.gmt" ) {
        if (file.exists( outfile) ) load( outfile )
        if (p$spatial.domain == "snowcrab" ) {
          id = bathymetry.db( DS="lookuptable.sse.snowcrab" )
          Z = Z[id,]
        }
        return (Z )
      }

      Z = bathymetry.db( p, DS="Z.planar" )
      dZ = bathymetry.db( p, DS="dZ.planar" )
      ddZ = bathymetry.db( p, DS="ddZ.planar" )
      Z = merge( Z, dZ, by=c("plon", "plat"), sort=FALSE )
      Z = merge( Z, ddZ, by=c("plon", "plat"), sort=FALSE )
      save (Z, file=outfile, compress=T )

			return( paste( "Completed:", outfile )  )
    }


  # ----------------


    if (DS %in% c("lookuptable.sse.snowcrab.redo", "lookuptable.sse.snowcrab" )) {
      #\\ DS="lookuptable.sse.snowcrab(.redo)" creates/returns a lookuptable for SSE -> snowcrab domains
      #\\   both share the same initial domains + resolutions and so it is faster to operate upon the indices
      fn = file.path( project.datadirectory("bio.bathymetry"), "interpolated", "sse.snowcrab.lookup.rdata")
      if (DS== "lookuptable.sse.snowcrab" ) {
        if (file.exists(fn)) load(fn)
        return(id)
      }
      zSSE = bathymetry.db ( p=spacetime_parameters( type="SSE" ), DS="baseline" )
      zSSE$id.sse = 1:nrow(zSSE)

      zsc  = bathymetry.db ( p=spacetime_parameters( type="snowcrab" ), DS="baseline" )
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

    if ( DS %in% c("bathymetry.spacetime.inputs.data", "bathymetry.spacetime.inputs.data.redo" )) {
      #\\ DS="bathymetry.spacetime.input" is a low-level call that prepares the bathymetry data
      #\\   for input into a table for further processing
      fn = file.path( datadir, paste( "bathymetry", "spacetime", "input", p$spatial.domain,  "rdata", sep=".") )
      if (DS =="bathymetry.spacetime.inputs.data" ) {
        load( fn)
        return( B )
      }
      print( "Warning: this needs a lot of RAM .. ~40GB depending upon resolution of discretization" )
      B = bathymetry.db ( p=p, DS="z.lonlat.rawdata" )
      B = B[ which(B$z > -100),]
      B = lonlat2planar( B, proj.type=p$internal.projection )
      B$plon = grid.internal( B$plon, p$plons )
      B$plat = grid.internal( B$plat, p$plats )
      B = block.spatial ( xyz=B[,c("plon", "plat", "z")], function.block=block.mean )
      save( B, file=fn, compress=TRUE)
      return(fn)
    }

    # --------------

    if ( DS %in% c("bathymetry.spacetime.inputs.prediction", "bathymetry.spacetime.inputs.prediction.redo" )) {
      #\\ DS="bathymetry.spacetime.input" is a low-level call that creates the input data table in a table
      fn = file.path( datadir, paste( "bathymetry", "spacetime", "input", "prediction", p$spatial.domain,  "rdata", sep=".") )
      if (DS =="bathymetry.spacetime.inputs.prediction" ) {
        load( fn)
        return( B )
      }
      B = expand.grid( p$plons, p$plats, KEEP.OUT.ATTRS=FALSE)
      names( B ) = c("plon", "plat")
      save( B, file=fn, compress=TRUE)
      return(fn)
    }


    # ----------------

    if ( DS == "landmasks.create" ) {
      # on resolution of predictions
      pps  =  expand.grid( plons=p$plons, plats=p$plats)
      V = SpatialPoints( planar2lonlat( pps, proj.type=p$internal.crs )[, c("lon", "lat" )], CRS("+proj=longlat +datum=WGS84") )
      landmask( lonlat=V, db="worldHires",regions=c("Canada", "US"), ylim=c(36,53), xlim=c(-72,-45), tag="predictions" )

      # on resolution of statistics
      V = as.data.frame(p$ff$Sloc[])
      names(V) = c("plon", "plat") 
      # V = data.frame( cbind(plon=Sloc[,1], plat=Sloc[,2]) )
      V = SpatialPoints( planar2lonlat( V, proj.type=p$internal.crs )[, c("lon", "lat" )], CRS("+proj=longlat +datum=WGS84") )
      landmask( lonlat=V, db="worldHires",regions=c("Canada", "US"), ylim=c(36,53), xlim=c(-72,-45), tag="statistics" )
    }

    #-------------------------

    if ( DS %in% c("bathymetry.spacetime.finalize.redo", "bathymetry.spacetime.finalize" )) {
      #// bathymetry( p, DS="bathymetry.spacetime.finalize(.redo)" return/create the
      #//   spacetime interpolated method formatted and finalised for production use
      fn = file.path(  project.datadirectory("bio.bathymetry"), "interpolated",
        paste( "bathymetry", "spacetime", "finalized", p$spatial.domain, "rdata", sep=".") )
      if (DS =="bathymetry.spacetime.finalize" ) {
        B = NULL
        if ( file.exists ( fn) ) load( fn)
        return( B )
      }

      nr = p$nplons
      nc = p$nplats

      BP = spacetime( p=p, DS="inla.predictions" )
      BP = BP[, c( "plon", "plat", "mean", "sdev")]
      names(BP) = c( "plon", "plat", "z", "Z.predictionSD") # really Z.mean but for historical compatibility "z"

      # remove land
      oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="land", tag="predictions" )
      BP$z[oc] = NA
      BP$Z.predictionSD[oc] = NA

      Bmn = matrix( data=BP$z, nrow=nr, ncol=nc )  # means

      # first order central differences but the central term drops out:
      # diffr = ( ( Bmn[ 1:(nr-2), ] - Bmn[ 2:(nr-1), ] ) + ( Bmn[ 2:(nr-1), ] - Bmn[ 3:nr, ] ) ) / 2
      # diffc = ( ( Bmn[ ,1:(nc-2) ] - Bmn[ ,2:(nc-1) ] ) + ( Bmn[ ,2:(nc-1) ] - Bmn[ ,3:nc ] ) ) / 2
      diffr =  Bmn[ 1:(nr-2), ] - Bmn[ 3:nr, ]
      diffc =  Bmn[ ,1:(nc-2) ] - Bmn[ ,3:nc ]
      rm (Bmn); gc()

      dZ = ( diffr[ ,2:(nc-1) ] + diffc[ 2:(nr-1), ] ) / 2
      dZ = rbind( dZ[1,], dZ, dZ[nrow(dZ)] )  # top and last rows are copies .. dummy value to keep dim correct
      dZ = cbind( dZ[,1], dZ, dZ[,ncol(dZ)] )

      BP$dZ =  abs(c(dZ))

      # gradients
      ddiffr =  dZ[ 1:(nr-2), ] - dZ[ 3:nr, ]
      ddiffc =  dZ[ ,1:(nc-2) ] - dZ[ ,3:nc ]
      rm( dZ ); gc()

      ddZ = ( ddiffr[ ,2:(nc-1) ] + ddiffc[ 2:(nr-1), ] ) / 2
      ddZ = rbind( ddZ[1,], ddZ, ddZ[nrow(ddZ)] )  # top and last rows are copies .. dummy value to keep dim correct
      ddZ = cbind( ddZ[,1], ddZ, ddZ[,ncol(ddZ)] )
      BP$ddZ = abs(c(ddZ))

      # merge into statistics
      BS = spacetime( p=p, DS="inla.statistics" )
      B = cbind( BP, BS )
      names(B) = c( names(BP), "Z.rangeMode", "Z.rangeSD", "Z.spatialSD", "Z.observationSD" )

      save( B, file=fn, compress=TRUE)
      return(fn)

      if (0) {
        aoi = which( B$z > 5 & B$z < 3000 & B$Z.rangeMode < 500)
        levelplot( log(z) ~ plon + plat, B[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( log(Z.rangeMode) ~ plon + plat, B[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( Z.rangeSD ~ plon + plat, B[aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )

      }
    }

    # ------------


    if ( DS %in% c( "complete", "complete.redo") ) {
      #// underlying storage data format is as a list of rasters
      #// regridding and selection to area of interest as specificied by girds.new=c("SSE", etc)
      Z = NULL

      if ( DS %in% c( "complete") ) {

        domain = NULL
        if ( is.null(domain)) {
          if ( exists("spatial.domain", p)) {
            domain = p$spatial.domain
          } else if ( exists( "grids.new", p) ) { # over-rides p$spatial domain
            if( length( p$grids.new )== 1 ) {
              domain = p$grids.new
        } } }

        fn = file.path( project.datadirectory("bio.bathymetry", "interpolated"),
          paste( "bathymetry", "spde_complete", domain, "rdata", sep=".") )
        if ( file.exists ( fn) ) load( fn)

        if ( return.format == "brick" ) {
          Z = brick(Z)
          return( Z )
        }

        if ( return.format == "dataframe" ) { ## default
          Z = as( brick(Z), "SpatialPointsDataFrame" )
          Z = as.data.frame(Z)
          u = names(Z)
          names(Z)[ which( u=="x") ] = "plon"
          names(Z)[ which( u=="y") ] = "plat"
          return( Z )
        }

        if ( return.format == "dataframe.filtered" ) {
          Z = as( brick(Z), "SpatialPointsDataFrame" )
          Z = as.data.frame(Z)
          u = names(Z)
          names(Z)[ which( u=="x") ] = "plon"
          names(Z)[ which( u=="y") ] = "plat"
          Z = filter.bathymetry( DS=p$spatial.domain, Z=Z )
          return( Z )
        }

        if ( return.format == "SpatialPointsDataFrame" ) { ## default
          Z = as( brick(Z), "SpatialPointsDataFrame" )
          return( Z )
        }
        if ( return.format %in% c("list") ) return( Z  )
      }

      p0 = p  # the originating parameters
      Z0 = bathymetry.db( p=p0, DS="bathymetry.spacetime.finalize" )
      coordinates( Z0 ) = ~ plon + plat
      crs(Z0) = crs( p0$interal.crs )
      above.sealevel = which( Z0$z < -1 ) # depth values < 0 are above  .. retain 1 m above to permits isobath calc
      if (length(above.sealevel)>0) Z0[ above.sealevel ] = NA

      Z = list()

      grids = unique( c( p$spatial.domain, p$grids.new ))

      for (gr in grids ) {
        p1 = spacetime_parameters( type=gr )
        for (vn in names(Z0)) {
          Z[[vn]] = projectRaster(
            from =rasterize( Z0, spacetime_parameters_to_raster(p0), field=vn, fun=mean),
            to   =spacetime_parameters_to_raster( p1) )
        }
        fn = file.path( project.datadirectory("bio.bathymetry", "interpolated"),
          paste( "bathymetry", "spde_complete", p1$spatial.domain, "rdata", sep=".") )
        save (Z, file=fn, compress=TRUE)
        print(fn)
      }
      return ( "Completed subsets" )
    }

  }  # end bathymetry.db








