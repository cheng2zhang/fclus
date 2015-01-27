/* Three-dimensional Lennard-Jones fluid */



"use strict";



/* initialize a fcc lattice */
function lj_initfcc3d(lj)
{
  var i, j, k, id, n = lj.n;

  var n1 = Math.floor(Math.pow(2*n, 1.0/3) + 0.999999); // # of particles per side
  var a = lj.l / n1;
  var noise = a * 1e-5;
  for (id = 0, i = 0; i < n1 && id < n; i++) {
    for (j = 0; j < n1 && id < n; j++) {
      for (k = 0; k < n1 && id < n; k++) {
        if ((i+j+k) % 2 === 0) {
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



/* get the tail correction */
function lj_gettail3d(lj, rho, n)
{
  var irc, irc3, irc6, utail, ptail;

  irc = 1 / lj.rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  utail = 8 * Math.PI * rho * n / 9 * (irc6 - 3) * irc3;
  ptail = 32 * Math.PI * rho * rho / 9 * (irc6 - 1.5) * irc3;
  return [utail, ptail];
}



/* annihilate the total angular momentum
 * solve
 *   /  y^2 + z^2    -x y      -x y      \
 *   |  -x y       X^2 + z^2   -y z      |  c  =  I
 *   \  -x z         -y z     x^2 + y^2  /
 * use a velocity field
 *    v = c X r
 *   */
function lj_shiftang3d(x, v, n)
{
  var i;
  var xc = [0, 0, 0], xi = [0, 0, 0], ang = [0, 0, 0], am = [0, 0, 0];
  var dv = [0, 0, 0];
  var mat = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  var xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;

  for (i = 0; i < n; i++) {
    vinc(xc, x[i]);
  }
  vsmul(xc, 1.0/n);
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(ang, xi, v[i]);
    vinc(am, ang);
    xx += xi[0]*xi[0];
    yy += xi[1]*xi[1];
    zz += xi[2]*xi[2];
    xy += xi[0]*xi[1];
    yz += xi[1]*xi[2];
    zx += xi[2]*xi[0];
  }
  mat[0][0] = yy+zz;
  mat[1][1] = xx+zz;
  mat[2][2] = xx+yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  var inv = rm3_inv(mat);
  ang[0] = -vdot(inv[0], am);
  ang[1] = -vdot(inv[1], am);
  ang[2] = -vdot(inv[2], am);
  // ang is the solution of M^(-1) * I
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross3d(dv, ang, xi);
    vinc(v[i], dv);
  }
}



var mousedown = 0;
var mousemoved = 0;
var mousex = -1;
var mousey = -1;
var viewmat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];



function ljmousedown(e)
{
  e = e || window.event;
  mousex = e.clientX;
  mousey = e.clientY;
  mousedown = 1;
  //console.log("mousedown", e.clientX, e.clientY, m2str(viewmat));
}



function ljmouseup(e)
{
  e = e || window.event;
  mousex = -1;
  mousey = -1;
  mousemoved = mousedown - 1;
  mousedown = 0;
  //console.log("mouseup", e.clientX, e.clientY, m2str(viewmat));
}



function ljmousemove(e)
{
  if ( !mousedown ) {
    return;
  }
  e = e || window.event;
  if ( mousex >= 0 && mousey >= 0 ) {
    var target = e.target ? e.target : e.srcElement;
    viewmat = mxrot3d(viewmat, 180.0 * (e.clientY - mousey) / target.height);
    viewmat = myrot3d(viewmat, 180.0 * (e.clientX - mousex) / target.width);
    paint(); // defined in ljmain.js
  }
  mousex = e.clientX;
  mousey = e.clientY;
  mousedown += 1;
  //console.log("mousemove", e.clientX, e.clientY, m2str(viewmat));
}



/* apply the view matrix */
function transform(x, l)
{
  var n = x.length;
  var xyz = newarr(n), xc = [l * 0.5, l * 0.5, l * 0.5], xi = [0, 0, 0];

  for ( var i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    xyz[i] = mmulv(viewmat, xi);
    vinc(xyz[i], xc);
  }
  return xyz;
}



function sortbyz(x)
{
  var i, j, k, l, n = x.length;
  var xyz = newarr2d(n, 3), rt = newarr(D);
  var idmap = newarr(n);
  var invmap = newarr(n);

  for ( i = 0; i < n; i++ ) {
    idmap[i] = i;
    // i:         index of the output array `xyz`
    // idmap[i]:  index of the input array `x`
    // so xyz[i] --> x[ idmap[i] ];
    invmap[i] = i;
  }

  // use bubble sort
  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[i]);
  }

  for ( i = 0; i < n; i++ ) {
    // find the ith smallest z
    k = i;
    var zmin = x[ idmap[i] ][2];
    for ( j = i + 1; j < n; j++ ) {
      if ( x[ idmap[j] ][2] < zmin ) {
        k = j;
        zmin = x[ idmap[j] ][2];
      }
    }
    if ( k != i ) {
      // before
      //  xyz[i] --> x[ idmap[i] ]
      //  xyz[k] --> x[ idmap[k] ]
      // after
      //  xyz[i] --> x[ idmap[k] ]
      //  xyz[k] --> x[ idmap[i] ]
      l = idmap[i];
      idmap[i] = idmap[k];
      idmap[k] = l;
    }
  }

  for ( i = 0; i < n; i++ ) {
    vcopy(xyz[i], x[ idmap[i] ]);
  }
  // compute the inverse map
  for ( i = 0; i < n; i++ ) {
    invmap[ idmap[i] ] = i;
  }
  return [xyz, idmap, invmap];
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



// return the z-scaling factor
function getzscale(r, zmin, zmax, ortho)
{
  if ( ortho ) {
    return 0.7;
  } else {
    var zf = (r[2] - zmin) / (zmax - zmin);
    return 0.7 + 0.3 * zf;
  }
}



// draw all atoms in the box
function ljdraw3d(lj, target, xin, userscale, edges, colors)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // the system dimension is L + two radii
  var scale = userscale * Math.min(width, height) / (lj.l + 1.0);
  var ortho = grab("orthographic").checked;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  var ret, xyz, idmap, invmap, zmin, zmax;
  var i, j, ic;
  var xi, yi, scli, xj, yj, sclj;

  // draw the box
  var l = lj.l;
  var xbox = [[0,0,0], [0,0,l], [0,l,0], [0,l,l], [l,0,0], [l,0,l], [l,l,0], [l,l,l]];
  var xboxt = transform(xbox, l); // apply the rotation matrix
  ret = sortbyz(xboxt); // sort the corners by the z index
  xyz = ret[0];
  idmap = ret[1];
  invmap = ret[2];
  zmax = xyz[7][2];
  zmin = xyz[0][2];
  //console.log(invmap, xbox, xboxt);
  ctx.lineWidth = 1;
  ctx.strokeStyle = '#c0c0c0';
  var boxedges = [[0,1], [2,3], [4,5], [6,7],
                  [0,2], [1,3], [4,6], [5,7],
                  [0,4], [1,5], [2,6], [3,7]];

  for ( ic = 0; ic < boxedges.length; ic++ ) {
    i = invmap[ boxedges[ic][0] ];
    scli = scale * getzscale(xyz[i], zmin, zmax, ortho);
    xi = Math.floor(  (xyz[i][0] - lj.l * 0.5) * scli + width  * 0.5 );
    yi = Math.floor( -(xyz[i][1] - lj.l * 0.5) * scli + height * 0.5 );

    j = invmap[ boxedges[ic][1] ];
    sclj = scale * getzscale(xyz[j], zmin, zmax, ortho);
    xj = Math.floor(  (xyz[j][0] - lj.l * 0.5) * sclj + width  * 0.5 );
    yj = Math.floor( -(xyz[j][1] - lj.l * 0.5) * sclj + height * 0.5 );

    drawLine(ctx, xi, yi, xj, yj);
    //console.log(xi, yi, xj, yj, i, j, boxedges[ic], xyz[i], xyz[j], scli, sclj);
  }

  // start to draw particles
  var xt = transform(xin, lj.l); // apply the rotation matrix
  ret = sortbyz(xt); // sort particles by the z order
  xyz = ret[0];
  idmap = ret[1];
  invmap = ret[2];
  // xyz[i]           --> xt[ idmap[i] ]
  // xyz[ invmap[i] ] --> xt[ i ]
  zmax = xyz[lj.n - 1][2];
  zmin = xyz[0][2];

  // draw lines that were used to group clusters
  if ( edges ) {
    ctx.lineWidth = 2;
    ctx.strokeStyle = '#808080';
    for ( ic = 0; ic < edges.length; ic++ ) {
      i = invmap[ edges[ic][0] ];
      scli = scale * getzscale(xyz[i], zmin, zmax, ortho);
      xi = Math.floor(  (xyz[i][0] - lj.l * 0.5) * scli + width  * 0.5 );
      yi = Math.floor( -(xyz[i][1] - lj.l * 0.5) * scli + height * 0.5 );

      j = invmap[ edges[ic][1] ];
      sclj = scale * getzscale(xyz[j], zmin, zmax, ortho);
      xj = Math.floor(  (xyz[j][0] - lj.l * 0.5) * sclj + width  * 0.5 );
      yj = Math.floor( -(xyz[j][1] - lj.l * 0.5) * sclj + height * 0.5 );
      drawLine(ctx, xi, yi, xj, yj, '#aaaaaa', 8);
      drawLine(ctx, xi, yi, xj, yj, '#bbbbbb', 4);
      drawLine(ctx, xi, yi, xj, yj, '#cccccc', 2);
    }
  }

  // determine which particles to circle around
  var ccnt = newarr(lj.g.nc);
  var mark = newarr(lj.n);
  for ( i = 0; i < lj.n; i++ ) {
    ic = lj.g.cid[i];
    mark[i] = ( ccnt[ic] < 1 );
    ccnt[ic] = 1;
  }

  // draw each particle
  for ( i = 0; i < lj.n; i++ ) {
    var z = xyz[i][2];
    var zf = (z - zmin) / (zmax - zmin);
    // make closer particles larger
    var zscl = getzscale(xyz[i], zmin, zmax, ortho);
    var scl = scale * zscl;
    var x = Math.floor(  (xyz[i][0] - lj.l * 0.5) * scl + width  * 0.5 );
    var y = Math.floor( -(xyz[i][1] - lj.l * 0.5) * scl + height * 0.5 );
    var rz = Math.floor( 0.5 * scl );
    ic = lj.g.cid[ idmap[i] ];
    var color = colors[ic];
    if ( ic === 0 ) {
      color = seedcolor;
    }
    if ( mark[ idmap[i] ] ) {
      // circle around the first particle of the cluster
      //drawBall(ctx, x, y, rz, color,     5); // outer outline
      //drawBall(ctx, x, y, rz, "#f0f0f0", 2); // inner outline

      // add some shade around the root vertices
      paintBall(ctx, x, y, rz * 1.3, transpColor(color, 0.0), color, 0, 0, rz * 0.8);
    }
    var spotcolor = "#e0e0e0";
    zscl = getzscale(xyz[i], zmin, zmax, false);
    var color2 = darkenColor(color, zscl);
    paintBall(ctx, x, y, rz, color2, spotcolor);
  }
}


