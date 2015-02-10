


"use strict";



/* hybrid MC */
CaGo.prototype.dohmc = function(hmc, wl)
{
  var acc;
  var rmsd = this.getRMSD(this.x);
  // hmc.farr[0] is the previous RMSD
  var dv = wl.getv( rmsd ) - wl.getv( hmc.fdat[0] );
  var iarr = null, farr = [ rmsd, this.epot ];

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
    this.epot = farr[1];
  }
  this.rmsd = farr[0];
  return acc;
};



/* compute the RMSD with a trial move */
CaGo.prototype.getRMSD2 = function(x, i, xi)
{
  for ( var j = 0; j < this.n; j++ ) {
    vcopy(go.x1[j], x[j]);
  }
  vcopy(go.x1[i], xi);
  return this.getRMSD(go.x1);
}



/* Metropolis algorithm with an RMSD bias */
CaGo.prototype.metrormsd = function(amp, bet)
{
  var i = Math.floor(this.n * rand01());
  var xi = newarr(D);
  for ( var d = 0; d < D; d++ ) {
    xi[d] = this.x[i][d] + amp * (rand01() * 2 - 1);
  }
  var du = this.depot(this.x, i, xi);

  var rmsd0 = this.rmsd;
  var rmsd1 = this.getRMSD2(this.x, i, xi);
  var dv = wl.getv( rmsd1 ) - wl.getv( rmsd0 );

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
    //if ( Math.abs(this.epot - this.force(this.x, this.x1)) > 1e-6 ) {
    //  stopsimul();
    //  throw new Error("e mismatch " + this.epot + " " + this.force(this.x, this.x1) + " " + du);
    //}
    this.rmsd = rmsd1;
    return 1;
  } else {
    return 0;
  }
};



