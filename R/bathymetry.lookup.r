
bathymetry.lookup = function( p, locs, vnames="z" ) {
  varnames = unique( c( "plon", "plat", vnames ) ) 
  domain = bathymetry.db(p=p, DS="baseline", varnames=varnames )
  domain_map = lbm::array_map( "xy->1", domain[,c("plon","plat")], gridparams=p$gridparams )
  locs_map = lbm::array_map( "xy->1", locs, gridparams=p$gridparams )
  locs_index = match( locs_map, domain_map )
  out = domain[locs_index, vnames]
  return(out)
}


