/* core functions of the LJ (Lennard-Jones) object */



"use strict";



function lj_initfcc(lj)
{
  if ( lj.dim == 2 ) {
    return lj_initfcc2d(lj);
  } else if ( lj.dim == 3 ) {
    return lj_initfcc3d(lj);
  }
}



function lj_gettail(lj, rho, n)
{
  if ( lj.dim == 2 ) {
    return lj_gettail2d(lj, rho, n);
  } else if ( lj.dim == 3 ) {
    return lj_gettail3d(lj, rho, n);
  }
}



function lj_setrho(lj, rho)
{
  lj.rho = rho;
  lj.vol = lj.n / rho;
  lj.l = Math.pow(lj.vol, 1.0 / lj.dim);
  lj.rc = Math.min( lj.l / 2, lj.rcdef );
  lj.rc2 = lj.rc * lj.rc;
  var irc = 1 / lj.rc;
  var irc2 = irc * irc;
  var irc6 = irc2 * irc2 * irc2;
  lj.epot_shift = 4 * irc6 * (irc6 - 1);
  var ret = lj_gettail(lj, rho, lj.n); // to be defined in lj2d.h or lj3d.h
  lj.epot_tail = ret[0];
  lj.p_tail = ret[1];
}



/* remove the center of mass motion */
function lj_rmcom(x, dim, n)
{
  var i;
  var rc = newarr(dim);

  for ( i = 0; i < n; i++ ) {
    vinc(rc, x[i]);
  }
  vsmul(rc, 1.0 / n);
  for ( i = 0; i < n; i++ ) {
    vdec(x[i], rc);
  }
}



function lj_shiftang(x, v, n)
{
  if ( D === 2 ) {
    lj_shiftang2d(x, v, n);
  } else if ( D === 3 ) {
    lj_shiftang3d(x, v, n);
  }
}



function LJ(n, dim, rho, rcdef, rcls)
{
  var i, d;

  this.n = n;
  this.dim = dim;
  this.dof = n * dim - dim * (dim + 1) / 2;
  this.rcdef = rcdef;
  this.x = newarr2d(n, dim);
  this.v = newarr2d(n, dim);
  this.f = newarr2d(n, dim);
  this.x2 = newarr2d(n, dim);
  this.r2ij = newarr2d(n, n);
  this.r2i = newarr(n);

  lj_setrho(this, rho);
  lj_initfcc(this);

  // initialize random velocities
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < dim; d++ ) {
      this.v[i][d] = randgaus();
    }
  }

  lj_rmcom(this.v, dim, n);
  lj_shiftang(this.x, this.v, n);

  this.epot = 0;
  this.eps = 0;
  this.ep0 = 0;
  this.vir = 0;
  this.ekin = 0;

  this.g = new Graph(n);
  this.g2 = new Graph(n);
  if ( rcls === undefined || rcls === null ) {
    rcls = 1.6;
  }
  this.rcls = rcls;
  this.vcls = newarr(n + 1);
  this.chist = newarr(n + 1);
  this.chist_cnt = 0;
}



/* OO version of lj_setrho() */
LJ.prototype.setrho = function(rho)
{
  lj_setrho(this, rho);
};



function lj_pbcdist2(dx, a, b, l, invl)
{
  return vsqr( vpbc(vdiff(dx, a, b), l, invl) );
}



/* build a graph */
function lj_mkgraph(lj, g, rm)
{
  var i, j, n = lj.n, rm2 = rm * rm;

  g.empty();
  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      if ( lj.r2ij[i][j] < rm2 ) {
        g.link(i, j);
      }
    }
  }
  g.clus();
}



LJ.prototype.mkgraph = function(g)
{
  lj_mkgraph(this, g, this.rcls);
}



/* build a graph with the distances from k computed from r2i */
function lj_mkgraph2(lj, g, k, rm)
{
  var i, j, n = lj.n;
  var r2, rm2 = rm * rm;

  g.empty();
  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      if ( i === k ) {
        r2 = lj.r2i[j];
      } else if ( j === k ) {
        r2 = lj.r2i[i];
      } else {
        r2 = lj.r2ij[i][j];
      }
      if ( r2 < rm2 ) {
        g.link(i, j);
      }
    }
  }
  g.clus();
}



LJ.prototype.mkgraph2 = function(g, k)
{
  lj_mkgraph2(this, g, k, this.rcls);
}



/* compute the cluster energy
 * call lj_mkgraph() first */
function lj_eclus(lj, g)
{
  return lj.vcls[ g.csize[ g.cid[0] ] ];
}



/* compute force and virial, return energy */
LJ.prototype.energy_low = function(x, r2ij)
{
  var dx = newarr(this.dim), dr2, dr6, ep, vir, rc2 = this.rc2;
  var l = this.l, invl = 1 / l;
  var i, j, npr = 0, n = this.n;

  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      r2ij[i][j] = dr2;
      r2ij[j][i] = dr2;
      if (dr2 < rc2) {
        dr2 = 1 / dr2;
        dr6 = dr2 * dr2 * dr2;
        vir += dr6 * (48 * dr6 - 24); // f.r
        ep += 4 * dr6 * (dr6 - 1);
        npr++;
      }
    }
  }
  return [ep + this.epot_tail, ep,
    ep - npr * this.epot_shift, // shifted energy
    vir];
};



LJ.prototype.energy = function()
{
  ret = this.energy_low(this.x, this.r2ij);
  this.epot = ret[0];
  this.ep0  = ret[1];
  this.eps  = ret[2];
  this.vir  = ret[3];
  return this.epot;
};



/* compute force and virial, return energy */
LJ.prototype.force_low = function(x, f, r2ij)
{
  var dx = newarr(this.dim), fi = newarr(this.dim);
  var dr2, dr6, fs, ep, vir, rc2 = this.rc2;
  var l = this.l, invl = 1/l;
  var i, j, npr = 0, n = this.n;

  for (i = 0; i < n; i++) {
    vzero(f[i]);
  }
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      r2ij[i][j] = dr2;
      r2ij[j][i] = dr2;
      if (dr2 < rc2) {
        dr2 = 1 / dr2;
        dr6 = dr2 * dr2 * dr2;
        fs = dr6 * (48 * dr6 - 24); // f.r
        vir += fs; // f.r
        fs *= dr2; // f.r / r^2
        vsinc(fi, dx, fs);
        vsinc(f[j], dx, -fs);
        ep += 4 * dr6 * (dr6 - 1);
        npr++;
      }
    }
    vinc(f[i], fi);
  }
  return [ep + this.epot_tail, ep,
    ep - npr * this.epot_shift, // shifted energy
    vir];
};



LJ.prototype.force = function()
{
  var ret = this.force_low(this.x, this.f, this.r2ij);
  this.epot = ret[0];
  this.ep0  = ret[1];
  this.eps  = ret[2];
  this.vir  = ret[3];
  return this.epot;
};



/* compute pressure */
LJ.prototype.calcp = function(tp)
{
  return (this.dof * tp + this.vir) / (this.dim * this.vol) + this.p_tail;
};



/* velocity-verlet */
LJ.prototype.vv = function(dt)
{
  var i, n = this.n;
  var dth = dt * 0.5, l = this.l;

  for (i = 0; i < n; i++) { // VV part 1
    vsinc(this.v[i], this.f[i], dth);
    vsinc(this.x[i], this.v[i], dt);
    vwrap(this.x[i], l);
  }
  this.force();
  for (i = 0; i < n; i++) { // VV part 2
    vsinc(this.v[i], this.f[i], dth);
  }
};



/* compute the kinetic energy */
function lj_ekin(v, n)
{
  var i;
  var ek = 0;
  for ( i = 0; i < n; i++ ) {
    ek += vsqr( v[i] );
  }
  return ek/2;
}



/* exact velocity rescaling thermostat */
function lj_vrescale_low(v, n, dof, tp, dt)
{
  var i;
  var c = (dt < 700) ? Math.exp(-dt) : 0;
  var ek1 = lj_ekin(v, n);
  var r = randgaus();
  var r2 = randchisqr(dof - 1);
  var ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp / 2 - ek1)
      + 2 * r * Math.sqrt(c * (1 - c) * ek1 * tp / 2);
  ek2 = Math.max(ek2, 0.0);
  var s = Math.sqrt(ek2 / ek1);
  for (i = 0; i < n; i++) {
    vsmul(v[i], s);
  }
  return ek2;
}



LJ.prototype.vrescale = function(tp, dt)
{
  return lj_vrescale_low(this.v, this.n, this.dof, tp, dt);
};



/* position Langevin barostat, with coordinates only
 * NOTE: the first parameter is the degree of freedom
 * the scaling is r = r*s
 * set cutoff to half of the box */
LJ.prototype.langp0 = function(dt, tp, pext, ensx)
{
  var pint, amp, s, dlnv;
  var i;

  pint = this.calcp(tp);
  amp = Math.sqrt(2 * dt);
  dlnv = ((pint - pext) * this.vol / tp + 1 - ensx) * dt + amp * randgaus();
  s = Math.exp( dlnv / this.dim );
  this.vol *= Math.exp( dlnv );
  this.setrho(this.n / this.vol);
  for ( i = 0; i < this.n; i++ ) {
    vsmul(this.x[i], s);
  }
  this.force();
};



/* displace a random particle i, return i */
LJ.prototype.randmv = function(xi, amp)
{
  var i = Math.floor(rand01() * this.n), d;
  for ( d = 0; d < this.dim; d++ ) {
    xi[d] = this.x[i][d] + (rand01() * 2 - 1) * amp;
  }
  return i;
};



/* compute pair energy */
function lj_pair(dr2, rc2)
{
  if (dr2 < rc2) {
    var invr2 = 1 / dr2;
    var invr6 = invr2 * invr2 * invr2;
    var vir = invr6 * (48 * invr6 - 24); // f.r
    var u  = 4 * invr6 * (invr6 - 1);
    return [true, u, vir];
  } else {
    return [false, 0.0, 0.0];
  }
}



/* return the energy change from displacing x[i] to xi */
LJ.prototype.depot = function(i, xi)
{
  var j, n = this.n;
  var l = this.l, invl = 1/l, rc2 = this.rc2, u, vir, ret;
  var dx = newarr(this.dim), r2;

  u = 0.0;
  vir = 0.0;
  for ( j = 0; j < n; j++ ) { // pair
    if ( j === i ) {
      continue;
    }
    r2 = ( i < j ) ? this.r2ij[i][j] : this.r2ij[j][i];
    ret = lj_pair(r2, rc2);
    if ( ret[0] ) {
      u -= ret[1];
      vir -= ret[2];
    }
    r2 = lj_pbcdist2(dx, xi, this.x[j], l, invl);
    ret = lj_pair(r2, rc2);
    if ( ret[0] ) {
      u += ret[1];
      vir += ret[2];
    }
    this.r2i[j] = r2;
  }
  return [u, vir];
};



/* commit a particle displacement */
LJ.prototype.commit = function(i, xi, du, dvir, ucls)
{
  var j;
  vwrap( vcopy(this.x[i], xi), this.l );
  this.ep0 += du;
  this.epot += du;
  this.vir += dvir;
  for ( j = 0; j < i; j++ ) {
    this.r2ij[j][i] = this.r2i[j];
  }
  for ( j = i + 1; j < n; j++ ) {
    this.r2ij[i][j] = this.r2i[j];
  }
  this.ecls = ucls;
  this.g.copy( this.g2 );
  // enable the following line to check if
  // this.r2ij and this.g are correct
  //lj_mkr2ij(this, this.x, this.r2ij, true);
};



/* Metropolis algorithm */
LJ.prototype.metro = function(amp, bet)
{
  var acc = 0;
  var xi = newarr(this.dim);

  var i = this.randmv(xi, amp);
  var ret = this.depot(i, xi);
  var du = ret[0];
  var dvir = ret[1];

  lj.mkgraph2(lj.g2, i);
  var ucls = lj_eclus(lj, lj.g2);
  var ducls = ucls - lj_eclus(lj, lj.g);

  if ( du + ducls < 0 ) {
    acc = 1;
  } else {
    var r = rand01();
    acc = ( r < Math.exp( -bet * (du + ducls) ) );
  }
  if ( acc ) {
    this.commit(i, xi, du, dvir, ucls);
    return 1;
  }
  return 0;
};



LJ.prototype.chist_clear = function()
{
  this.chist_cnt = 0;
  for ( var i = 0; i < this.n; i++ ) {
    this.chist[i] = 0;
  }
};



LJ.prototype.chist_add = function(g)
{
  this.chist_cnt += 1;
  this.chist[ g.csize[ g.cid[0] ] ] += 1;
};



LJ.prototype.chist_tostr = function()
{
  var s = "";

  if ( this.chist_cnt > 0 ) {
    for ( var i = 0; i <= this.n; i++ ) {
      if ( this.chist[i] > 0 ) {
        s += "" + i + " " + (this.chist[i]/this.chist_cnt) + " " + this.vcls[i] + "\n";
      }
    }
  }
  return s;
};



LJ.prototype.clus = function()
{
  this.mkgraph(this.g);
  this.chist_add(this.g);
};



LJ.prototype.update_vcls = function(g, lnf)
{
  var i;
  this.vcls[ g.csize[ g.cid[0] ] ] += lnf;

  // find the minimal of the potential and subtract it
  var min = 1e30;
  for ( i = 1; i <= this.n; i++ ) {
    min = Math.min(this.vcls[i], min);
  }
  for ( i = 1; i <= this.n; i++ ) {
    this.vcls[i] -= min;
  }
};



/* check the flatness of the histogram, return the new lnf */
LJ.prototype.update_lnf = function(lnf, flatness, frac)
{
  var i, n = this.n;
  var sh = 0.0, shh = 0.0, h;

  for ( i = 1; i <= n; i++ ) {
    h = this.chist[i];
    sh += h;
    shh += h * h;
  }
  if ( sh <= 0 ) return lnf;
  sh /= n;
  shh = Math.max(shh / n - sh * sh, 0);

  this.hflatness = Math.sqrt(shh) / sh;
  if ( this.hflatness < flatness ) {
    this.chist_clear();
    console.log("change lnf from " + lnf + " to " + (lnf * frac) );
    lnf *= frac;
  }
  return lnf;
};



/* update r2ij */
function lj_mkr2ij(lj, x, r2ij, check)
{
  var dx = newarr(lj.dim), dr2, rc2 = lj.rc2, rm2 = lj.rcls * lj.rcls;
  var l = lj.l, invl = 1 / l;
  var i, j, n = lj.n;

  for ( i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if ( check ) {
        if ( Math.abs( r2ij[i][j] - dr2 ) > 1e-5 ) {
          console.log("r2ij corrupted for i", i, " j ", j, ": ", dr2, " ", r2ij[i][j]);
          throw new Error("r2ij corruption");
          stop();
        }
        var cij = (dr2 < rm2);
        if ( lj.g.linked(i, j) != cij ) {
          console.log("graph corrupted for i", i, " j ", j, ": ", dr2, " ", r2ij[i][j], " rm2 " , rm2, " cij ", cij, " gij ", lj.g.linked(i, j));
          throw new Error("cij corruption");
          stop();
        }
      }
      r2ij[i][j] = dr2;
      r2ij[j][i] = dr2;
    }
  }
}



/* wrap coordinates to group particles in the same cluster */
function lj_wrapclus(lj, xin, xout, g, rm)
{
  var ic, i, j, n = lj.n, head, end;
  var l = lj.l, invl = 1 / l, dx = newarr(lj.dim);
  var edges = [];

  for ( i = 0; i < n; i++ ) {
    vcopy(xout[i], xin[i]);
  }

  // shift xout[0] to the center of the box
  vdiff(dx, [l * 0.5, l * 0.5, l * 0.5], xout[0]);
  for ( i = 0; i < n; i++ ) {
    vinc( xout[i], dx );
    vwrap( xout[i], l );
  }

  //lj_mkr2ij(lj, xout, lj.r2ij, true);
  lj.mkgraph(g);

  for ( ic = 0; ic < g.nc; ic++ ) {
    // find the seed of this cluster
    for ( i = 0; i < n; i++ ) {
      if ( g.cid[i] === ic ) {
        break;
      }
    }
    if ( i >= n ) {
      throw new Error("no particle in cluster " + ic);
      return null;
    }

    g.queue[ head = 0 ] = i;
    g.cid[ i ] = -1;
    end = 1;
    for ( ; head < end; head++ ) {
      i = g.queue[ head ];
      // add neighbors of i into the queue
      for ( j = 0; j < n; j++ ) {
        if ( g.linked(i, j) && g.cid[j] >= 0 ) {
          if ( g.cid[j] !== ic ) {
            throw new Error("cluster id of " + j + " " + g.cid[j] + " " + ic);
          }
          g.cid[j] = -1;
          g.queue[ end++ ] = j;
          vpbc( vdiff(dx, xout[j], xout[i]), l, invl );
          vinc( vcopy(xout[j], xout[i]), dx );
          edges.push( [i, j] );
        }
      }
    }
    if ( end != g.csize[ic] ) {
      throw new Error("cluster " + ic + ": size " + end + " vs " + g.csize[ic]);
    }
  }

  return edges;
}




/* list edges */
function lj_listedges(lj, x)
{
  var i, j, n = lj.n;
  var edges = [], dx = newarr(lj.dim), dr2, rm2 = lj.rcls * lj.rcls;

  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      dr2 = vsqr( vdiff(dx, x[i], x[j]) );
      // if i and j are properly wrapped
      // we add a bond between them
      if ( dr2 < rm2 ) {
        edges.push( [i, j] );
      }
    }
  }
  return edges;
}

