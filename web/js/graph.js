/* Graph */



"use strict";




function Graph(n)
{
  this.n = n;
  this.mat = newarr2d(n, n);
  this.nc = 0;
  this.csize = newarr(n);
  this.cseed = newarr(n);
  this.cid = newarr(n);
  this.queue = newarr(n);
  this.icmax = -1;
};



Graph.prototype.empty = function()
{
  var n = this.n;
  for ( var i = 0; i < n; i++ ) {
    for ( var j = 0; j < n; j++ ) {
      this.mat[i][j] = false;
    }
  }
}



Graph.prototype.link = function(i, j)
{
  this.mat[i][j] = true;
  this.mat[j][i] = true;
}



Graph.prototype.unlink = function(i, j)
{
  this.mat[i][j] = false;
  this.mat[j][i] = false;
}



Graph.prototype.linked = function(i, j)
{
  return this.mat[i][j];
}



/* g = g2 */
Graph.prototype.copy = function(g2)
{
  var i, j, ic, n = g2.n;

  if ( this.n != n ) {
    throw new Error("cannot copy graphs of different sizes " + this.n + " vs " + n);
    return -1;
  }
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < n; j++ )
      this.mat[i][j] = g2.mat[i][j];

  // copy the cluster information
  this.nc = g2.nc;
  for ( i = 0; i < n; i++ ) {
    this.csize[i] = g2.csize[i];
    this.cid[i] = g2.cid[i];
  }
  for ( ic = 0; ic < this.nc; ic++ ) {
    this.cseed[ic] = g2.cseed[ic];
  }
  return 0;
}



/* do clustering */
Graph.prototype.clus = function(seed)
{
  var i, j, n = this.n, ic;
  var head, end, max;

  for ( i = 0; i < n; i++ ) this.cid[i] = -1;

  max = 0;
  for ( ic = 0; ; ic++ ) {
    if ( ic === 0 && seed !== undefined && seed !== null ) {
      // seed of the first cluster
      i = seed;
    } else {
      // find the first free atom, whose cid is negative
      for ( i = 0; i < n; i++ )
        if ( this.cid[i] < 0 )
          break;
      if ( i >= n ) break;
    }

    // grow a cluster starting from i
    this.queue[ head = 0 ] = i;
    this.cid[i] = ic;
    this.cseed[ic] = i;
    end = 1;
    for ( ; head < end; head++ ) {
      i = this.queue[head];
      // add neighbors of i into the queue
      for ( j = 0; j < n; j++ ) {
        if ( this.linked(i, j) && this.cid[j] < 0 ) {
          this.cid[j] = ic;
          this.queue[ end++ ] = j;
        }
      }
    }
    // now `end` is the number of particles in cluster `ic`
    this.csize[ic] = end;
    if ( end > max ) {
      max = end;
      this.icmax = ic;
    }
  }
  this.nc = ic;
  for ( ic = 0; ic < this.nc; ic++ ) {
    i = this.cseed[ic];
    var ic2 = this.cid[i];
    if ( ic2 !== ic ) {
      console.log("ic", ic, "ic2", ic2, "seed", i);
    }
  }
}



Graph.prototype.clus_print = function()
{
  var s = "" + this.nc + " clusters: ";
  for ( var i = 0; i < this.nc; i++ ) {
    s += "" + this.csize[i] + " ";
  }
  return s;
}



function checkgraph(g)
{
  console.log(g.nc, "clusters");
  for ( var ic = 0; ic < g.nc; ic++ ) {
    var i = lj.g.cseed[ic];
    var ic2 = lj.g.cid[i];
    if ( ic != ic2 ) {
      console.log("ic", ic, "seed", i, "cid", ic2);
      throw new Error("bad cluster");
    }
  }
}



