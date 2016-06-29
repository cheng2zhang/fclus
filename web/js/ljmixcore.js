/* core functions of the LJ (Lennard-Jones) object */



"use strict";



/* initialize a fcc lattice */
function ljmix_initfcc2d(lj)
{
  var i, j, id = 0, n = lj.n;
  var n1 = Math.floor( Math.sqrt(2*n) + 0.999999 ); // # of particles per side
  var a = lj.l / n1;
  var noise = a * 1e-5;

  for ( id = 0, i = 0; i < n1 && id < n; i++ ) {
    for ( j = 0; j < n1 && id < n; j++ ) {
      if ( (i + j) % 2 === 0 ) {
        // add some noise to prevent two atoms happened to
        // be separated precisely by the cutoff distance,
        // which might be half of the box
        lj.x[id][0] = (i + 0.5) * a + noise * (2*rand01() - 1);
        lj.x[id][1] = (j + 0.5) * a + noise * (2*rand01() - 1);
        id++;
      }
    }
  }
}



/* initialize a fcc lattice */
function ljmix_initfcc3d(lj)
{
  var i, j, k, id, n = lj.n;

  var n1 = Math.floor(Math.pow(2*n, 1.0/3) + 0.999999); // # of particles per side
  var a = lj.l / n1;
  var noise = a * 1e-5;
  for (id = 0, i = 0; i < n1 && id < n; i++) {
    for (j = 0; j < n1 && id < n; j++) {
      for (k = 0; k < n1 && id < n; k++) {
        if ((i + j + k) % 2 === 0) {
          /* add some noise to prevent two atoms happened to
           * be separated by precisely some special cutoff distance,
           * which might be half of the box */
          lj.x[id][0] = (i + 0.5) * a + noise * (2*rand01() - 1);
          lj.x[id][1] = (j + 0.5) * a + noise * (2*rand01() - 1);
          lj.x[id][2] = (k + 0.5) * a + noise * (2*rand01() - 1);
          id++;
        }
      }
    }
  }
}



function ljmix_initfcc(lj)
{
  if ( lj.dim == 2 ) {
    return ljmix_initfcc2d(lj);
  } else if ( lj.dim == 3 ) {
    return ljmix_initfcc3d(lj);
  }
}



/* get the tail correction */
function ljmix_gettail2d(lj, rho, n)
{
  var is, js, ns = lj.ns;
  var rd, irc, irc3, irc6, sig, eps, utail, ptail;

  rd = lj.rc * lj.rc;
  utail = 0;
  ptail = 0;
  for ( is = 0; is < ns; is++ ) {
    for ( js = 0; js < ns; js++ ) {
      sig = lj.sigij[is][js];
      eps = lj.epsij[is][js];
      irc = sig / lj.rc;
      irc3 = irc * irc * irc;
      irc6 = irc3 * irc3;
      utail += Math.PI*eps*lj.rho[is]*lj.rho[js]*(0.4*irc6 - 1)*irc6*rd;
      ptail += Math.PI*eps*lj.rho[is]*lj.rho[js]*(2.4*irc6 - 3)*irc6*rd;
    }
  }
  return [utail, ptail];
}



/* get the tail correction */
function ljmix_gettail3d(lj, rho, n)
{
  var is, js, ns = lj.ns;
  var rd, irc, irc3, irc6, sig, eps, utail, ptail;

  rd = lj.rc * lj.rc * lj.rc;
  utail = 0;
  ptail = 0;
  for ( is = 0; is < ns; is++ ) {
    for ( js = 0; js < ns; js++ ) {
      sig = lj.sigij[is][js];
      eps = lj.epsij[is][js];
      irc = sig / lj.rc;
      irc3 = irc * irc * irc;
      irc6 = irc3 * irc3;
      utail +=  8*Math.PI*eps*lj.rho[is]*lj.rho[js]/9*(irc6 - 3.0)*irc6*rd;
      ptail += 32*Math.PI*eps*lj.rho[is]*lj.rho[js]/9*(irc6 - 1.5)*irc6*rd;
    }
  }
  return [utail, ptail];
}



/* apply the view matrix */

function ljmix_gettail(lj, rho, n)
{
  if ( lj.dim == 2 ) {
    return ljmix_gettail2d(lj, rho, n);
  } else if ( lj.dim == 3 ) {
    return ljmix_gettail3d(lj, rho, n);
  }
}



function ljmix_setvol(lj, vol)
{
  var is;

  lj.vol = vol;
  for ( is = 0; is < lj.ns; is++ ) {
    lj.rho[is] = lj.np[is] / vol;
  }
  lj.l = Math.pow(lj.vol, 1.0 / lj.dim);
  lj.rc = Math.min( lj.l * 0.5, lj.rcdef );
  lj.rc2 = lj.rc * lj.rc;
  var ret = ljmix_gettail(lj, rho, lj.n);
  lj.epot_tail = ret[0];
  lj.p_tail = ret[1];
}



function LJMix(np, sig, eps, dim, rho, rcdef, rcls, vcls)
{
  var i, j, t, n, is, js, d;
  var vol;

  this.ns = np.length;
  this.np = np;
  this.sig = sig;
  this.eps = eps;

  // compute the total number of particles
  n = 0;
  for ( is = 0; is < this.ns; is++ ) {
    n += this.np[is];
  }
  this.n = n;

  // compute the densities
  vol = this.n / rho;
  this.rho = newarr( this.ns );
  for ( is = 0; is < this.ns; is++ ) {
    this.rho[is] = this.np[is] / vol;
  }

  // compute the pairwise LJ parameters
  this.sigij = newarr2d(this.ns, this.ns);
  this.epsij = newarr2d(this.ns, this.ns);
  for ( is = 0; is < this.ns; is++ ) {
    for ( js = 0; js < this.ns; js++ ) {
      this.sigij[is][js] = (this.sig[is] + this.sig[js]) / 2;
      this.epsij[is][js] = Math.sqrt( this.eps[is] * this.eps[js] );
    }
  }

  // assign the types
  this.type = newarr(this.n);
  i = 0;
  for ( is = 0; is < this.ns; is++ ) {
    for ( j = 0; j < this.np[is]; j++ ) {
      this.type[ i++ ] = is;
    }
  }

  this.dim = dim;
  this.dof = n * dim - dim;
  this.rcdef = rcdef;
  this.x = newarr2d(n, dim);
  this.v = newarr2d(n, dim);
  this.f = newarr2d(n, dim);
  this.x2 = newarr2d(n, dim);
  this.r2ij = newarr2d(n, n);
  this.r2i = newarr(n);

  ljmix_setvol(this, vol);

  ljmix_initfcc(this);

  // randomly swap coordinates to mix things up
  for ( t = 0; t < n * n; t++ ) {
    i = Math.floor(rand01() * n);
    j = (i + 1 + Math.floor(rand01() * (n - 1))) % n;
    vswap(this.x[i], this.x[j]);
  }

  // initialize random velocities
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < dim; d++ ) {
      this.v[i][d] = randgaus();
    }
  }

  rmcom(this.v, null, n);
  shiftang(this.x, this.v, null, n);

  this.epot = 0;
  this.ep0 = 0;
  this.vir = 0;
  this.ekin = 0;

  var np0 = this.np[0];

  this.g = new Graph(np0);
  this.g2 = new Graph(np0);
  if ( rcls === undefined || rcls === null ) {
    rcls = 1.6;
  }
  this.rcls = rcls;
  this.vcls = vcls;
  this.chistall = newarr(np0);
  this.chistall_cnt = 0;
  this.cseed = 0; // seed of the cluster

  this.renv = newarr(n); // envelope radii
}



/* OO version of ljmix_setvol() */
LJMix.prototype.setvol = function(rho)
{
  ljmix_setvol(this, rho);
};



function ljmix_pbcdist2(dx, a, b, l, invl)
{
  return vsqr( vpbc(vdiff(dx, a, b), l, invl) );
}



/* build a graph */
function ljmix_mkgraph(lj, g, rm)
{
  var i, j, np = lj.np[0], rm2 = rm * rm;

  g.empty();
  for ( i = 0; i < np; i++ ) {
    for ( j = i + 1; j < np; j++ ) {
      if ( lj.r2ij[i][j] < rm2 ) {
        g.link(i, j);
      }
    }
  }
  g.clus(lj.cseed);
}



LJMix.prototype.mkgraph = function(g)
{
  ljmix_mkgraph(this, g, this.rcls);
}



/* build a graph with the distances from k computed from r2i */
function ljmix_mkgraph2(lj, g, k, rm)
{
  var i, j, np = lj.np[0];
  var r2, rm2 = rm * rm;

  g.empty();
  for ( i = 0; i < np; i++ ) {
    for ( j = i + 1; j < np; j++ ) {
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
  g.clus(lj.cseed);
}



LJMix.prototype.mkgraph2 = function(g, k)
{
  ljmix_mkgraph2(this, g, k, this.rcls);
}



/* compute the cluster energy
 * call mkgraph() first */
function ljmix_eclus(lj, g)
{
  return lj.vcls[ g.csize[ g.cid[ lj.cseed ] ] - 1 ];
}



/* compute force and virial, return energy */
LJMix.prototype.energy_low = function(x, r2ij)
{
  var dx = newarr(this.dim), dr2, dr6, ep, vir, rc2 = this.rc2;
  var l = this.l, invl = 1 / l, sig, eps;
  var i, j, itp, jtp, npr = 0, n = this.n;

  ep = vir = 0;
  for ( i = 0; i < n - 1; i++ ) {
    itp = this.type[i];
    for ( j = i + 1; j < n; j++ ) {
      jtp = this.type[j];
      dr2 = ljmix_pbcdist2(dx, x[i], x[j], l, invl);
      r2ij[i][j] = dr2;
      r2ij[j][i] = dr2;
      if ( dr2 >= rc2 ) {
        continue;
      }

      sig = this.sigij[itp][jtp];
      eps = this.epsij[itp][jtp];
      dr2 = (sig * sig) / dr2;
      dr6 = dr2 * dr2 * dr2;
      vir += eps * dr6 * (48 * dr6 - 24); // f.r
      ep += eps * 4 * dr6 * (dr6 - 1);
      npr++;
    }
  }
  return [ep + this.epot_tail, ep, vir];
};



LJMix.prototype.energy = function()
{
  ret = this.energy_low(this.x, this.r2ij);
  this.epot = ret[0];
  this.ep0  = ret[1];
  this.vir  = ret[2];
  return this.epot;
};



/* compute force and virial, return energy */
LJMix.prototype.force_low = function(x, f, r2ij)
{
  var dx = newarr(this.dim), fi = newarr(this.dim);
  var dr2, dr6, fs, ep, vir, rc2 = this.rc2;
  var l = this.l, invl = 1/l, sig, eps;
  var i, j, itp, jtp, npr = 0, n = this.n;

  for (i = 0; i < n; i++) {
    vzero(f[i]);
  }
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    itp = this.type[i];
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = ljmix_pbcdist2(dx, x[i], x[j], l, invl);
      r2ij[i][j] = dr2;
      r2ij[j][i] = dr2;
      if ( dr2 >= rc2 ) {
        continue;
      }
      jtp = this.type[j];
      sig = this.sigij[itp][jtp];
      eps = this.epsij[itp][jtp];
      dr2 = (sig * sig) / dr2;
      dr6 = dr2 * dr2 * dr2;
      fs = eps * dr6 * (48 * dr6 - 24); // f.r
      vir += fs; // f.r
      fs *= dr2; // f.r / r^2
      vsinc(fi, dx, fs);
      vsinc(f[j], dx, -fs);
      ep += eps * 4 * dr6 * (dr6 - 1);
      npr++;
    }
    vinc(f[i], fi);
  }
  return [ep + this.epot_tail, ep, vir];
};



LJMix.prototype.force = function()
{
  var ret = this.force_low(this.x, this.f, this.r2ij);
  this.epot = ret[0];
  this.ep0  = ret[1];
  this.vir  = ret[2];
  return this.epot;
};



/* compute pressure */
LJMix.prototype.calcp = function(tp)
{
  return (this.dof * tp + this.vir) / (this.dim * this.vol) + this.p_tail;
};



/* velocity-verlet */
LJMix.prototype.vv = function(dt)
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



/* exact velocity rescaling thermostat */
LJMix.prototype.vrescale = function(tp, dt)
{
  return md_vrescale(this.v, null, this.n, this.dof, tp, dt);
};



/* position Langevin barostat, with coordinates only
 * NOTE: the first parameter is the degree of freedom
 * the scaling is r = r*s
 * set cutoff to half of the box */
LJMix.prototype.langp0 = function(dt, tp, pext, ensx)
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
LJMix.prototype.randmv = function(xi, amp)
{
  var i = Math.floor(rand01() * this.n), d;
  for ( d = 0; d < this.dim; d++ ) {
    xi[d] = this.x[i][d] + (rand01() * 2 - 1) * amp;
  }
  return i;
};



/* compute pair energy */
function ljmix_pair(dr2, sig, eps, rc2)
{
  if (dr2 < rc2) {
    var invr2 = (sig * sig) / dr2;
    var invr6 = invr2 * invr2 * invr2;
    var vir = eps * invr6 * (48 * invr6 - 24); // f.r
    var u  = eps * 4 * invr6 * (invr6 - 1);
    return [true, u, vir];
  } else {
    return [false, 0.0, 0.0];
  }
}



/* return the energy change from displacing x[i] to xi */
LJMix.prototype.depot = function(i, xi)
{
  var j, itp, jtp, n = this.n;
  var l = this.l, invl = 1/l, rc2 = this.rc2, u, vir, ret;
  var dx = newarr(this.dim), r2, sig, eps;

  u = 0.0;
  vir = 0.0;
  itp = this.type[i];
  for ( j = 0; j < n; j++ ) { // pair
    if ( j === i ) {
      continue;
    }
    r2 = ( i < j ) ? this.r2ij[i][j] : this.r2ij[j][i];
    jtp = this.type[j];
    sig = this.sigij[ itp ][ jtp ];
    eps = this.epsij[ itp ][ jtp ];
    ret = ljmix_pair(r2, sig, eps, rc2);
    if ( ret[0] ) {
      u -= ret[1];
      vir -= ret[2];
    }
    r2 = ljmix_pbcdist2(dx, xi, this.x[j], l, invl);
    ret = ljmix_pair(r2, sig, eps, rc2);
    if ( ret[0] ) {
      u += ret[1];
      vir += ret[2];
    }
    this.r2i[j] = r2;
  }
  return [u, vir];
};



/* commit a particle displacement */
LJMix.prototype.commit = function(i, xi, du, dvir, ucls)
{
  var j;
  vwrap( vcopy(this.x[i], xi), this.l );
  this.ep0 += du;
  this.epot += du;
  this.vir += dvir;
  for ( j = 0; j < i; j++ ) {
    this.r2ij[j][i] = this.r2i[j];
  }
  for ( j = i + 1; j < this.n; j++ ) {
    this.r2ij[i][j] = this.r2i[j];
  }
  this.ecls = ucls;
  this.g.copy( this.g2 );
  // enable the following line to check if
  // this.r2ij and this.g are correct
  //ljmix_mkr2ij(this, this.x, this.r2ij, true);
};



/* Metropolis algorithm
 * graph this.g is updated */
LJMix.prototype.metro = function(amp, bet)
{
  var acc = 0;
  var xi = newarr(this.dim);

  var i = this.randmv(xi, amp);
  var ret = this.depot(i, xi);
  var du = ret[0];
  var dvir = ret[1];

  this.mkgraph2(this.g2, i);
  var ucls = ljmix_eclus(lj, this.g2);
  var ducls = ucls - ljmix_eclus(this, this.g);
  var dutot = bet * du + ducls;

  if ( dutot < 0 ) {
    acc = 1;
  } else {
    var r = rand01();
    acc = ( r < Math.exp( -dutot ) );
  }
  if ( acc ) {
    this.commit(i, xi, du, dvir, ucls);
    return 1;
  }
  return 0;
};



/* change the seed for clustering */
LJMix.prototype.changeseed = function(g)
{
  var np = this.np[0];
  var i = (this.cseed + 1 + Math.floor(rand01() * (np - 1))) % np;
  var sz0 = g.csize[ g.cid[ this.cseed ] ];
  var sz1 = g.csize[ g.cid[ i ] ];
  var acc;
  if ( sz0 === sz1 ) {
    acc = true;
  } else {
    var dv = this.vcls[ sz1 - 1 ] - this.vcls[ sz0 - 1 ];
    if ( dv < 0 ) {
      acc = true;
    } else {
      var r = rand01();
      acc = ( r < Math.exp(-dv) );
    }
  }
  if ( acc ) {
    this.cseed = i;
    this.mkgraph(this.g);
  }
};



/* hybrid MC */
LJMix.prototype.dohmc = function(hmc)
{
  var acc;
  // compute the current cluster size
  var csize = this.g.csize[ this.g.cid[this.cseed] ];
  // hmc.iarr[0] is the previous size
  var dv = this.vcls[ csize - 1 ] - this.vcls[ hmc.idat[0] - 1 ];
  var iarr = [ csize, this.cseed ];
  var farr = [ this.epot ];

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
    this.cseed = iarr[1];
    this.epot = farr[0];
    //this.force();
    //this.mkgraph(this.g);
  }
  this.csize = iarr[0];
  return acc;
};



/* update r2ij */
function ljmix_mkr2ij(lj, x, r2ij, check)
{
  var dx = newarr(lj.dim), dr2, rc2 = lj.rc2, rm2 = lj.rcls * lj.rcls;
  var l = lj.l, invl = 1 / l;
  var i, j, n = lj.n;

  for ( i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = ljmix_pbcdist2(dx, x[i], x[j], l, invl);
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
function ljmix_wrapclus(lj, xin, xout, g)
{
  var ic, i, j, n = lj.n, head, end;
  var l = lj.l, invl = 1 / l, dx = newarr(lj.dim);
  var np = lj.np[0];
  var mat = newarr2d(np, np);

  for ( i = 0; i < n; i++ ) {
    vcopy(xout[i], xin[i]);
  }

  // shift xout[lj.cseed] to the center of the box
  vdiff(dx, [l * 0.5, l * 0.5, l * 0.5], xout[lj.cseed]);
  //console.log("before", lj.cseed, xout[lj.cseed], dx);
  for ( i = 0; i < n; i++ ) {
    vinc( xout[i], dx );
    vwrap( xout[i], l );
  }

  //ljmix_mkr2ij(lj, xout, lj.r2ij, true);
  lj.mkgraph(g);
  //console.log("after", lj.cseed, xout[lj.cseed], g.cid[lj.cseed], lj.g.cseed[0]);

  var gg = lj.g;
  if ( gg.cid[lj.cseed] !== 0 || gg.cseed[0] !== lj.cseed ) {
    stop();
    throw new Error("bad seed! " + lj.cseed + " " + gg.cseed[0] + " " + g.cseed[0]);
  }

  for ( ic = 0; ic < g.nc; ic++ ) {
    // find the seed of this cluster
    i = g.cseed[ic];

    g.queue[ head = 0 ] = i;
    g.cid[ i ] = -1;
    end = 1;
    for ( ; head < end; head++ ) {
      i = g.queue[ head ];
      // add neighbors of i into the queue
      for ( j = 0; j < np; j++ ) {
        if ( g.linked(i, j) && g.cid[j] >= 0 ) {
          if ( g.cid[j] !== ic ) {
            throw new Error("cluster id of " + j + " " + g.cid[j] + " " + ic);
          }
          g.cid[j] = -1;
          g.queue[ end++ ] = j;
          vpbc( vdiff(dx, xout[j], xout[i]), l, invl );
          vinc( vcopy(xout[j], xout[i]), dx );
          mat[i][j] = 1;
          mat[j][i] = 1;
        }
      }
    }
    if ( end != g.csize[ic] ) {
      throw new Error("cluster " + ic + ": size " + end + " vs " + g.csize[ic]);
    }
  }

  return mat;
}




/* get the clustering adjacency matrix */
function ljmix_getclsmat(lj, x)
{
  var i, j, np = lj.np[0];
  var dx = newarr(lj.dim), dr2, rm2 = lj.rcls * lj.rcls;

  var mat = newarr2d(np, np);
  for ( i = 0; i < np; i++ ) {
    for ( j = i + 1; j < np; j++ ) {
      dr2 = vsqr( vdiff(dx, x[i], x[j]) );
      // if i and j are properly wrapped
      // we add a bond between them
      if ( dr2 < rm2 ) {
        mat[i][j] = 1;
        mat[j][i] = 1;
      }
    }
  }
  return mat;
}




