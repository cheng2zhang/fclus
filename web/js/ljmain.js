/* Handle web interface */



"use strict";



var lj = null;
var wl = null;
var hmc = null;

var n = 55;
var rho = 0.7;
var tp = 1.5;
var rcdef = 1000.0;
var rcls = 1.6;

var timer_interval = 100; // in milliseconds
var ljtimer = null;
var simulmethod = "MD";

var mddt = 0.002;
var thdt = 0.02;
var nstepspsmd = 100; // number of steps per second for MD
var nstepspfmd = 10;  // number of steps per frame for MD
var nstepsmd = 0;
var hmctot = 0.0;
var hmcacc = 0.0;
var nsthmc = 1;
var nvswaps = 1;

var mcamp = 0.2;
var nstepspsmc = 10000; // number of steps per second for MC
var nstepspfmc = 1000;  // number of steps per frame for MC
var nstepsmc = 0;
var mctot = 0.0;
var mcacc = 0.0;

var wl_lnf0 = 0.01;
var wl_flatness = 0.3;
var wl_frac = 0.5;
var invt_c = 1.0;

var changeseed = true;

var histplot = null;
var vplot = null;

var adjustscale = 1.0;
var xpaint = null; // coordinates used for painting
var randcolors = null; // random colors for each cluster
var seedcolor = "#ff2010"; // color of the special cluster



function getparams()
{
  n = get_int("n", 55);
  var dim = get_int("dimension", 2);
  if ( dim === 2 || dim === 3 ) {
    D = dim;
  }
  rho = get_float("density", 0.7);
  tp = get_float("temperature", 1.5);
  rcdef = get_float("rcutoff", 1000.0);
  rcls = get_float("rcluster", 1.6);

  timer_interval = get_int("timer_interval");
  simulmethod = grab("simulmethod").value;
  mddt = get_float("mddt", 0.002);
  thdt = get_float("thermostatdt", 0.01);
  nsthmc = get_int("nsthmc", 1);
  nvswaps = get_int("nvswaps", 1);
  nstepspsmd = get_int("nstepspersecmd", 100);
  nstepspfmd = nstepspsmd * timer_interval / 1000;
  nstepsmd = 0;

  mcamp = get_float("mcamp", 0.2);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;
  nstepsmc = 0;

  wl_lnf0 = get_float("wl_lnfinit", 0.01);
  wl_flatness = get_float("wl_flatness", 0.3);
  wl_frac = get_float("wl_frac", 0.5);
  invt_c = get_float("invt_c", 1.0);

  changeseed = grab("changeseed").checked;

  mousescale = get_float("ljscale");

  randcolors = newarr(n + 1);
  for ( var ic = 0; ic <= n; ic++ ) {
    randcolors[ic] = randHueColor(40, 120);
  }
}



function changescale()
{
  mousescale = get_float("ljscale");
  paint();
}



/* return a string of the current simulation details */
function getsinfo()
{
  var s = "", clist = "";

  var flatness = wl.getflatness();
  s += '<span class="math">ln <i>f</i> </span>: ' + wl.lnf.toExponential(3) + ".<br>";
  s += 'flatness: ' + roundto(flatness * 100, 2) + "%.<br>";
  s += 'seed: ' + lj.cseed + ".<br>";
  s += 'clusters: ';
  for ( var ic = 0; ic < lj.g.nc; ic++ ) {
    clist += "" + lj.g.csize[ic];
    if ( ic < lj.g.nc - 1 ) {
      clist += ", ";
    } else {
      clist += ".";
    }
  }
  var nlen = clist.length;
  while ( nlen < 25 ) {
    clist += "&nbsp;";
    nlen += 1;
  }
  s += clist + "<br>";
  return s;
}



function domd()
{
  var istep, sinfo = "";

  for ( istep = 1; istep <= nstepspfmd; istep++ ) {
    lj.vv(mddt);
    lj.vrescale(tp, thdt);
    lj.mkgraph(lj.g);
    if ( changeseed ) {
      lj.changeseed(lj.g);
    }

    // use hybrid MC to sample a flat histogram
    if ( nsthmc > 0 && istep % nsthmc === 0 ) {
      hmctot += 1;
      hmcacc += lj.dohmc(hmc);
      lj.ekin = md_vscramble(lj.v, null, lj.n, nvswaps);
    } else {
      lj.csize = lj.g.csize[ lj.g.cid[ lj.cseed ] ];
    }

    wl.add( lj.csize );
    wl.updatelnf();
  }
  wl.trimv();
  nstepsmd += nstepspfmd;
  sinfo += 'step ' + nstepsmd + ", ";
  sinfo += "hmcacc: " + roundto(100.0 * hmcacc / hmctot, 2) + "%.<br>";
  sinfo += getsinfo();
  return sinfo;
}



function domc()
{
  var istep, sinfo = "";

  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    // do a step of Metropolis algorithm
    mctot += 1.0;
    mcacc += lj.metro(mcamp, 1.0 / tp);
    // try to change the seed of the cluster, unnecessary
    if ( changeseed ) {
      lj.changeseed(lj.g);
    }
    lj.csize = lj.g.csize[ lj.g.cid[ lj.cseed ] ];
    wl.add( lj.csize );
    if ( wl.updatelnf() ) {
      // adjust the MC move size, to make acceptance ratio close to 0.5
      if ( grab("adjmcamp").checked ) {
        var avacc = mcacc/mctot;
        var acctarget = 0.5
        if ( Math.abs(avacc - acctarget) > 0.02 ) {
          mcamp *= Math.max( Math.min( Math.sqrt(avacc/acctarget), 2 ), 0.5 );
        }
        grab("mcamp").value = roundto(mcamp, 4);
      }
      console.log("acc", mcacc/mctot, ", amp", mcamp);
      mctot = 1e-30;
      mcacc = 0;
    }
  }
  wl.trimv();
  nstepsmc += nstepspfmc;
  sinfo += "step: " + nstepsmc + ", ";
  sinfo += "acc: " + roundto(100.0 * mcacc / mctot, 2) + "%.<br>";
  sinfo += getsinfo();
  return sinfo;
}



// normalize the histogram such that the maximum is 1.0
function normalize_hist(arr, n)
{
  var hs = newarr(n);
  var i, m = 0;
  for ( i = 0; i < n; i++ ) {
    m = Math.max(m, 1.0 * arr[i]);
  }
  for ( i = 0; i < n; i++ ) {
    hs[i] = 1.0 * arr[i] / m;
  }
  return hs;
}



/* update the histogram plot */
function updatehistplot(wl)
{
  var i;
  var dat = "Cluster size,Histogram (all time),Histogram (this stage)\n";

  var chhs = normalize_hist(wl.hh, wl.n);
  var chs  = normalize_hist(wl.h, wl.n);
  for ( i = 0; i < lj.n; i++ )
    dat += "" + (i+1) + "," + chhs[i] + "," + chs[i] + "\n";
  if ( histplot === null ) {
    var h = grab("ljbox").height / 2 - 5;
    var w = h * 3 / 2;
    var options = {
      //title: 'Histogram of cluster size',
      xlabel: '<small>Cluster size, <i>s</i></small>',
      ylabel: '<small>Histogram, <i>H</i></small>',
      includeZero: true,
      drawPoints: true,
      axisLabelFontSize: 10,
      pointSize: 2,
      xRangePad: 2,
      plotter: barChartPlotter,
      width: w,
      height: h
    };
    histplot = new Dygraph(grab("histplot"), dat, options);
  } else {
    histplot.updateOptions({ file: dat });
  }
}



/* update the cluster potential plot */
function updatevplot(wl)
{
  var i;
  var dat = "Cluster size,potential\n";
  for ( i = 0; i < wl.n; i++ )
    dat += "" + (i+1) + "," + wl.v[i] + "\n";
  if ( vplot === null ) {
    var h = grab("ljbox").height / 2 - 5;
    var w = h * 3 / 2;
    var options = {
      //title: 'Adaptive potential',
      xlabel: '<small>Cluster size, <i>s</i></small>',
      ylabel: '<small>Adaptive potential, <i>V</i></small>',
      includeZero: true,
      drawPoints: true,
      axisLabelFontSize: 10,
      pointSize: 2,
      xRangePad: 2,
      width: w,
      height: h
    };
    vplot = new Dygraph(grab("vplot"), dat, options);
  } else {
    vplot.updateOptions({ file: dat });
  }
}



function paint()
{
  if ( !lj ) {
    return;
  }

  var groupclus = grab("groupcluster").checked;
  var paintmat = null; // adjacency matrix for visualization
  if ( groupclus ) {
    xpaint = lj.x2;
    paintmat = lj_wrapclus(lj, lj.x, xpaint, lj.g2);
    // if we group particles according to clusters
    // some particles will flow out of the box
    // so we need a smaller scale
    adjustscale = 0.7;
  } else {
    xpaint = lj.x;
    paintmat = lj_getclsmat(lj, lj.x);
    adjustscale = 1.0;
  }

  var s = mousescale * adjustscale;
  var ballscale = get_float("ballscale", 1.0);

  if ( lj.dim === 2 ) {
    ljdraw2d(lj, "ljbox", xpaint, s, paintmat, randcolors, ballscale);
  } else if ( lj.dim === 3 ) {
    ljdraw3d(lj, "ljbox", xpaint, s, paintmat, randcolors, ballscale);
  }
}



function pulse()
{
  var sinfo;

  // the following parameter might have been changed
  if ( wl ) {
    wl.flatness = get_float("wl_flatness", 0.3);
    wl.frac = get_float("wl_frac", 0.5);
    wl.c = get_float("invt_c", 1.0);
  }

  if ( simulmethod === "MD" ) {
    sinfo = domd();
  } else if ( simulmethod === "MC" ) {
    sinfo = domc();
  }
  lj_clusvol(lj, lj.g);
  sinfo += "cluster volume: " + roundto(lj.clsvol[0], 2) + "/" + roundto(lj.vol, 2) + "\n";
  grab("sinfo").innerHTML = sinfo;

  paint();

  updatehistplot(wl);
  updatevplot(wl);
}



function stopsimul()
{
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
  }
  grab("pause").innerHTML = "&#9724;";
  hmctot = 1e-30;
  hmcacc = 0.0;
  mctot = 1e-30;
  mcacc = 0.0;
  munit(viewmat);
}



function pausesimul()
{
  if ( !lj ) return;
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
    grab("pause").innerHTML = "&#10704";
  } else {
    ljtimer = setInterval(
        function() { pulse(); },
        timer_interval);
    grab("pause").innerHTML = "&#9724;";
  }
}



function startsimul()
{
  stopsimul();
  getparams();
  wl = new WL(1, n + 1, 1, false, wl_lnf0, wl_flatness, wl_frac, invt_c, 0);
  lj = new LJ(n, D, rho, rcdef, rcls, wl.v);
  lj.force();
  lj.mkgraph(lj.g);
  hmc = new HMC(lj.n, 1, 1);
  hmc.push(lj.x, lj.v, lj.f,
      [ lj.g.csize[ lj.g.cid[ lj.cseed ] ] ], [lj.epot]);
  installmouse("ljbox", "ljscale");
  ljtimer = setInterval(
    function(){ pulse(); },
    timer_interval);
}



// on a mouse click
function pausesimul2()
{
  // skip a mouse-move
  if ( mousemoved > 0 ) {
    return;
  }
  if ( !lj ) {
    startsimul();
  } else if ( mousemoved === 0 ) {
    pausesimul();
  }
}



function changeparams()
{
  if ( ljtimer !== null ) {
    startsimul();
  }
}



function showtab(who)
{
  who = grab(who);
  var par = who.parentNode;
  var c = par.childNodes;
  var i, iwho, k = 0;

  // arrange the tabs
  for ( i = 0; i < c.length; i++ ) {
    if ( c[i].className === "params-panel" ) {
      if ( c[i] !== who ) {
        c[i].style.zIndex = k;
        k += 1;
      } else {
        iwho = k;
      }
    }
  }
  who.style.zIndex = k;

  // arrange the clickable tab titles
  k += 1;
  var pt = grab("tabsrow");
  pt.style.zIndex = k;
  var ct = pt.childNodes, ik = 0;
  for ( i = 0; i < ct.length; i++ ) {
    if ( ct[i].tagName ) {
      if ( ik === iwho ) {
        ct[i].style.fontWeight = "bold";
        ct[i].style.borderTop = "2px solid #c0c0d0";
      } else {
        ct[i].style.fontWeight = "normal";
        ct[i].style.borderTop = "0px solid #e0e0f0";
      }
      ik += 1;
    }
  }
}



function resizecontainer(a)
{
  var canvas = grab("ljbox");
  var ctx = canvas.getContext("2d");
  var w, h;
  if ( a === null || a === undefined ) {
    w = canvas.width;
    h = canvas.height;
  } else {
    a = parseInt( grab(a).value );
    w = h = a;
    canvas.width = w;
    canvas.height = h;
  }
  ctx.font = "24px Verdana";
  ctx.fillText("Click to start", w/2-40, h/2-10);

  var hsbar = 30; // height of the global scaling bar
  var hcbar = 40; // height of the control bar
  var htbar = 30; // height of the tabs bar
  var wr = h*3/4; // width of the plots
  var wtab = w; // width of the tabs
  var htab = 240;

  grab("simulbox").style.width = "" + w + "px";
  grab("simulbox").style.height = "" + h + "px";
  grab("simulbox").style.top = "" + hsbar + "px";
  grab("controlbox").style.top = "" + (h + hsbar) + "px";
  grab("ljscale").style.width = "" + (w - 100) + "px";
  histplot = null;
  grab("histplot").style.left = "" + w + "px";
  grab("histplot").style.width = "" + wr + "px";
  grab("vplot").style.top = "" + hcbar + "px";
  grab("histplot").style.height = "" + h/2 + "px";
  vplot = null;
  grab("vplot").style.left = "" + w + "px";
  grab("vplot").style.width = "" + wr + "px";
  grab("vplot").style.top = "" + (h/2 + hcbar) + "px";
  grab("vplot").style.height = "" + h/2 + "px";
  grab("tabsrow").style.top = "" + (h + hsbar + hcbar) + "px";
  grab("tabsrow").style.width = "" + wtab + "px";
  var c = grab("container").childNodes;
  var i;
  /* tabs */
  for ( i = 0; i < c.length; i++ ) {
    if ( c[i].className === "params-panel" ) {
      c[i].style.top = "" + (h + hsbar + hcbar + htbar) + "px";
      c[i].style.width = "" + (w - 20) + "px";
      c[i].style.height = "" + htab + "px";
    }
  }
  grab("sinfo").style.top = "" + (h + hsbar + hcbar + htbar) + "px";
  grab("sinfo").style.left = "" + (w + 10) + "px";
  grab("container").style.height = "" + (h + hsbar + hcbar + htbar + htab) + "px";
  grab("container").style.width = "" + (w + wr) + "px";
}



function init()
{
  resizecontainer();
  showtab("system-params");
}

