/* hybrid MC */



"use strict";



function HMC(n, ni, nf)
{
  this.n = n;
  this.x = newarr(n);
  this.v = newarr(n);
  this.f = newarr(n);
  this.ni = ni;
  if ( ni > 0 ) {
    this.idat = newarr(ni);
  }
  this.nf = nf;
  if ( nf > 0 ) {
    this.fdat = newarr(nf);
  }
}



HMC.prototype.push = function(x, v, f, idat, fdat)
{
  var n = this.n;
  cparr(this.x, x, n);
  cparr(this.v, v, n);
  cparr(this.f, f, n);
  if ( idat ) {
    cparr(this.idat, idat, this.ni);
  }
  if ( fdat ) {
    cparr(this.fdat, fdat, this.nf);
  }
};



HMC.prototype.pop = function(x, v, f, idat, fdat, reversev)
{
  var n = this.n;

  if ( reversev ) {
    for ( var i = 0; i < n; i++ ) {
      vneg( this.v[i] );
    }
  }
  cparr(x, this.x, n);
  cparr(v, this.v, n);
  cparr(f, this.f, n);
  if ( idat ) {
    cparr(idat, this.idat, this.ni);
  }
  if ( fdat ) {
    cparr(fdat, this.fdat, this.nf);
  }
}



