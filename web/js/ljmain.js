/* Handle web interface */



"use strict";



var lj = null;
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

var mcamp = 0.2;
var nstepspsmc = 10000; // number of steps per second for MC
var nstepspfmc = 1000;  // number of steps per frame for MC
var mctot = 0.0;
var mcacc = 0.0;

var sum1 = 1e-30;
var sumU = 0.0;
var sumP = 0.0;

var wl_lnf = 0.01;
var wl_flatness = 0.3;
var wl_frac = 0.5;

var histplot = null;
var vclsplot = null;

var userscale = 1.0;
var xpaint = null; // the coordinates to paint
var paintedges = null;
var randcolors = null;



function getparams()
{
  n = get_int("n", 55);
  var dim = get_int("dimension", 2);
  if ( dim === 2 || dim === 3 ) {
    D = dim;
  }
  rho = get_float("density", 0.7);
  temp = get_float("temp", 1.5);
  rcdef = get_float("rcutoff", 1000.0);
  rcls = get_float("rcluster", 1.6);

  simulmethod = grab("simulmethod").value;
  mddt = get_float("mddt", 0.002);
  thdt = get_float("thermostatdt", 0.01);
  nstepspsmd = get_int("nstepspersecmd", 100);
  nstepspfmd = nstepspsmd * timer_interval / 1000;

  mcamp = get_float("mcamp", 0.2);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;

  wl_lnf = get_float("wl_lnfinit", 0.01);
  wl_flatness = get_float("wl_flatness", 0.3);
  wl_frac = get_float("wl_frac", 0.5);

  randcolors = newarr(n + 1);
  for ( var ic = 0; ic <= n; ic++ ) {
    randcolors[ic] = randhuecolor(40, 120);
  }
  // set the first color to red
  randcolors[0] = "#ff2010";
}



/* for mouse wheel event */
function ljwheel(e){
  var delta = 0; // positive for scrolling up
  e = e || window.event;
  if ( e.wheelDelta ) { // IE/Opera
    delta = e.wheelDelta / 120;
  } else if ( e.detail ) { // Firefox
    delta = -e.detail / 3;
  }
  if ( delta > 0 ) {
    userscale *= 1.05;
  } else if ( delta < 0 ) {
    userscale *= 0.95;
  }
  //console.log("wheel", delta);
  if ( e.preventDefault ) {
    e.preventDefault();
  }
  e.returnValue = false;
  paint(); // defined later
}



/* for the wheel event */
function installwheel(target, handler)
{
  if ( target.addEventListener ) {
    // for IE9+, Chrome, Safari, Opera
    target.addEventListener('mousewheel', handler, false);
    // for Firefox
    target.addEventListener('DOMMouseScroll', handler, false);
  } else { // for IE 6/7/8
    target.attachEvent("onmousewheel", handler);
  }
}



function installmouse()
{
  var target = grab("ljbox");
  target.onmousedown = ljmousedown;
  target.onmouseup = ljmouseup;
  target.onmousemove = ljmousemove;
  installwheel(target, ljwheel);
}



function domd()
{
  var istep, sinfo = "";

  for ( istep = 0; istep < nstepspfmd; istep++ ) {
    lj.vv(mddt);
    lj.vrescale(tp, thdt);
    sum1 += 1.0;
    sumU += lj.epot / lj.n;
    sumP += lj.calcp(tp);
  }
  sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: ' + roundto(sumU/sum1, 3) + ", ";
  sinfo += '<span class="math"><i>P</i></span>: ' + roundto(sumP/sum1, 3);
  return sinfo;
}



function domc()
{
  var istep, sinfo = "";

  for ( istep = 0; istep < nstepspfmc; istep++ ) {
    mctot += 1.0;
    mcacc += lj.metro(mcamp, 1.0 / tp);
    lj.chist_add(lj.g);
    lj.update_vcls(lj.g, wl_lnf);
    wl_lnf = lj.update_lnf(wl_lnf, wl_flatness, wl_frac);
    sum1 += 1.0;
    sumU += lj.epot / lj.n;
    sumP += lj.calcp(tp);
  }
  sinfo += "acc: " + roundto(100.0 * mcacc / mctot, 2) + "%, ";
  sinfo += '<span class="math">ln <i>f</i> </span>: ' + wl_lnf + ", ";
  sinfo += 'flatness: ' + roundto(lj.hflatness * 100, 2) + "%. ";
  sinfo += '<br>clusters: ';
  for ( var ic = 0; ic < lj.g.nc; ic++ ) {
    sinfo += "" + lj.g.csize[ic];
    if ( ic < lj.g.nc - 1 ) {
      sinfo += ", ";
    } else {
      sinfo += ".";
    }
  }
  //sinfo += '<span class="math"><i>U</i>/<i>N</i></span>: ' + roundto(sumU/sum1, 3) + ", ";
  //sinfo += '<span class="math"><i>P</i></span>: ' + roundto(sumP/sum1, 3);
  return sinfo;
}



/* update the histogram plot */
function updatehistplot(lj)
{
  var i;
  var dat = "Cluster size,Histogram\n";
  for ( i = 1; i <= lj.n; i++ )
    dat += "" + i + "," + (lj.chist[i] / lj.chist_cnt) + "\n";
  if ( histplot === null ) {
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
      width: 360,
      height: 240
    };
    histplot = new Dygraph(grab("histplot"), dat, options);
  } else {
    histplot.updateOptions({ file: dat });
  }
}



/* update the cluster potential plot */
function updatevclsplot(lj)
{
  var i;
  var dat = "Cluster size,potential\n";
  for ( i = 1; i <= lj.n; i++ )
    dat += "" + i + "," + lj.vcls[i] + "\n";
  if ( vclsplot === null ) {
    var options = {
      //title: 'Adaptive potential',
      xlabel: '<small>Cluster size, <i>s</i></small>',
      ylabel: '<small>Adaptive potential, <i>V</i></small>',
      includeZero: true,
      drawPoints: true,
      axisLabelFontSize: 10,
      pointSize: 2,
      xRangePad: 2,
      width: 360,
      height: 240
    };
    vclsplot = new Dygraph(grab("vclsplot"), dat, options);
  } else {
    vclsplot.updateOptions({ file: dat });
  }
}



function paint()
{
  if ( !lj ) return;
  if ( lj.dim === 2 ) {
    ljdraw2d(lj, "ljbox", xpaint, userscale, paintedges, randcolors);
  } else if ( lj.dim === 3 ) {
    ljdraw3d(lj, "ljbox", xpaint, userscale, paintedges, randcolors);
  }
}



function pulse()
{
  var sinfo;

  if ( simulmethod === "MD" ) {
    sinfo = domd();
  } else if ( simulmethod === "MC" ) {
    sinfo = domc();
  }
  grab("sinfo").innerHTML = sinfo;

  var groupclus = grab("groupcluster").checked;
  paintedges = null;
  if ( groupclus ) {
    xpaint = lj.x2;
    paintedges = lj_wrapclus(lj, lj.x, xpaint, lj.g2, lj.rcls);
  } else {
    xpaint = lj.x;
    paintedges = lj_listedges(lj, lj.x); // list edges
  }

  paint();

  updatehistplot(lj);
  updatevclsplot(lj);
  //console.log( lj.chist_tostr() );
}



function stopsimul()
{
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
  }
  lj = null;
  mctot = 0.0;
  mcacc = 0.0;
  sum1 = 1e-30;
  sumU = 0.0;
  sumP = 0.0;
  munit(viewmat);
}



function pausesimul()
{
  if ( !lj ) return;
  if ( ljtimer !== null ) {
    clearInterval(ljtimer);
    ljtimer = null;
    grab("pause").value = "Resume";
  } else {
    ljtimer = setInterval(
        function() { pulse(); },
        timer_interval);
    grab("pause").value = "Pause";
  }
}



function startsimul()
{
  stopsimul();
  getparams();
  lj = new LJ(n, D, rho, rcdef, rcls);
  lj.force();
  lj.mkgraph(lj.g);
  installmouse();
  ljtimer = setInterval(
    function(){ pulse(); },
    timer_interval);
}



function changeparams()
{
  if ( ljtimer !== null ) {
    startsimul();
  }
}

