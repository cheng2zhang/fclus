


"use strict";



/* hybrid MC */
CaGo.prototype.dohmcnc = function(hmc, wl)
{
  var acc;
  var nc = this.ncontacts(this.x);
  // hmc.iarr[0] is the previous nc
  var dv = wl.getv( nc ) - wl.getv( hmc.idat[0] );
  var iarr = [ nc ], farr = [ this.epot ];

  if ( dv <= 0 ) {
    acc = 1;
  } else {
    var r = rand01();
    acc = ( r < Math.exp( -dv ) );
  }
  if ( acc ) {
    hmc.push(this.x, this.v, this.f, iarr, farr);
  } else {
    hmc.pop(this.x, this.v, this.f, iarr, farr, true);
    this.epot = farr[0];
  }
  this.nc = iarr[0];
  return acc;
};



/* compute the number of contacts with a trial move */
CaGo.prototype.ncontacts2 = function(x, i, xi)
{
  for ( var j = 0; j < this.n; j++ ) {
    vcopy(this.x1[j], x[j]);
  }
  vcopy(this.x1[i], xi);
  return this.ncontacts(this.x1);
}



/* Metropolis algorithm with an RMSD bias */
CaGo.prototype.metronc = function(amp, bet)
{
  var i = Math.floor(this.n * rand01());
  var xi = newarr(D);
  for ( var d = 0; d < D; d++ ) {
    xi[d] = this.x[i][d] + amp * (rand01() * 2 - 1);
  }
  var du = this.depot(this.x, i, xi);

  var nc0 = this.nc;
  var nc1 = this.ncontacts2(this.x, i, xi);
  var dv = wl.getv( nc1 ) - wl.getv( nc0 );

  var dutot = bet * du + dv;
  var acc;

  if ( dutot < 0 ) {
    acc = 1;
  } else {
    var r = rand01();
    acc = ( r < Math.exp( -dutot ) );
  }
  if ( acc ) {
    vcopy(this.x[i], xi);
    this.epot += du;
    this.nc = nc1;
    return 1;
  } else {
    return 0;
  }
};



