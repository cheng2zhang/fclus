/* miscellaneous routines */



"use strict";



// compute the mutual volume of two spheres of radii ri and rj
// with the distance of dij
function mutvol2d(ri, rj, dij)
{
  var ai = Math.acos( (ri*ri + dij*dij - rj*rj) / (2*ri*dij) );
  var vi = ri * ri * (ai - Math.sin(ai));
  var aj = Math.acos( (rj*rj + dij*dij - ri*ri) / (2*rj*dij) );
  var vj = rj * rj * (aj - Math.sin(aj));
  return vi + vj;
}



// compute the mutual volume of two spheres of radii ri and rj
// with the mutual distance of dij
function mutvol3d(ri, rj, dij)
{
  var xi = (ri*ri + dij*dij - rj*rj)/(2*dij);
  var vi = Math.PI * ( 2 * ri * ri * ri / 3 - ri * ri * xi + xi * xi * xi / 3 );
  var xj = (rj*rj + dij*dij - ri*ri)/(2*dij);
  var vj = Math.PI * ( 2 * rj * rj * rj / 3 - rj * rj * xj + xj * xj * xj / 3 );
  return vi + vj;
}



/* find the envelope radius */
function lj_getrenv(lj, g)
{
  var i, j, n = lj.n, ic, jc;
  var dr2, dr2min;

  // NOTE: currently, we'll use the uniform radius for all balls
  // for simplicity.
  for ( i = 0; i < n; i++ ) {
    lj.renv[i] = lj.rcls * 0.5;
  }
  return lj.renv;

  for ( i = 0; i < n; i++ ) {
    // find the closest particle that is not in this cluster
    ic = g.cid[i];
    dr2min = 1e30;
    for ( j = 0; j < n; j++ ) {
      jc = g.cid[j];
      if ( jc === ic ) continue;
      if ( i < j ) {
        dr2 = lj.r2ij[i][j];
      } else {
        dr2 = lj.r2ij[j][i];
      }
      if ( dr2 < dr2min ) {
        dr2min = dr2;
      }
    }
    lj.renv[i] = Math.sqrt( dr2min ) * 0.5;
  }

  // TODO: refine these parameters here
  var zinmax = (D == 3) ? 6 : 3;
  // exclude internal points
  for ( i = 0; i < n; i++ ) {
    ic = g.cid[i];
    var r2i = lj.renv[i] * lj.renv[i];
    var zin = 0;
    for ( j = 0; j < n; j++ ) {
      jc = g.cid[j];
      if ( jc !== ic || j === i ) continue;
      if ( i < j ) {
        dr2 = lj.r2ij[i][j];
      } else {
        dr2 = lj.r2ij[j][i];
      }
      if ( dr2 < r2i ) {
        zin += 1;
      }
    }
    // reduce radius for the internal points
    if ( zin > zinmax ) {
      lj.renv[i] = rcls * .5;
    }
  }
  //console.log(lj.renv);
}



// compute the cluster volume
function lj_clusvol(lj, g)
{
  var ic, i, j, n = lj.n;
  var vol = 0, ri, rj, dij;

  lj_getrenv(lj, g);
  lj.clsvol = newarr(g.nc);

  for ( ic = 0; ic < g.nc; ic++ ) {
    for ( i = 0; i < n; i++ ) {
      if ( g.cid[i] !== ic ) {
        continue;
      }
      ri = lj.renv[i];
      if ( D === 2 ) {
        vol += Math.PI * ri * ri;
      } else if ( D === 3 ) {
        vol += 4 * Math.PI / 3 * ri * ri * ri;
      }
      // compute the mutual exclusion volume
      for ( j = i + 1; j < n; j++ ) {
        var rj = lj.renv[j];
        var dij = Math.sqrt( lj.r2ij[i][j] );
        if ( dij >= ri + rj
          || ri + dij >= rj
          || rj + dij >= ri ) continue;
        if ( D === 2 ) {
          vol -= mutvol2d(ri, rj, dij);
        } else {
          vol -= mutvol3d(ri, rj, dij);
        }
        if ( isNaN(vol) ) {
          console.log( ic, i, j, ri, rj, dij, vol );
          throw new Error("invalid volume");
        }
      }
    }
    lj.clsvol[ic] = vol;
  }

  return lj.clsvol;
}




