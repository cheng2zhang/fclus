/* graphics routines */



"use strict";



/* draw line */
function drawLine(ctx, xi, yi, xj, yj)
{
  ctx.beginPath();
  ctx.moveTo(xi, yi);
  ctx.lineTo(xj, yj);
  ctx.stroke();
}



/* draw a ball that is centered at (x, y) with radius r
 * color is the color of the ball
 * the format of color is "#rrggbb" */
function drawBall(ctx, x, y, r, color, lineWidth)
{
  if ( color  ) {
    ctx.strokeStyle = color;
  }
  if ( lineWidth ) {
    ctx.lineWidth = lineWidth;
  }
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2*Math.PI);
  ctx.closePath();
  ctx.stroke();
}



/* draw a ball that is centered at (x, y) with radius r
 * color is the color of the ball
 * spotcolor is the color of the spotlight
 * the format of color is "#rrggbb" */
function paintBall(ctx, x, y, r, color, spotcolor,
    spotx, spoty, spotr)
{
  if ( spotcolor === undefined || spotcolor === null ) {
    spotcolor = "#a0a0a0";
  }
  if ( spotx === undefined || spotcolor === null ) {
    spotx = r * 0.3;
  }
  if ( spoty === undefined || spoty === null ) {
    spoty = r * 0.4;
  }
  if ( spotr === undefined || spotr === null ) {
    spotr = r * 0.1;
  }
  var grd = ctx.createRadialGradient(x + spotx, y - spoty, spotr, x, y, r);
  grd.addColorStop(0, spotcolor); // spotlight color
  grd.addColorStop(1, color); // ball color
  ctx.fillStyle = grd;
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2*Math.PI);
  ctx.closePath();
  ctx.fill();
}



function rgb2str(r, g, b)
{
  r = Math.floor(r).toString(16);
  if ( r.length == 1 ) {
    r = "0" + r;
  }

  g = Math.floor(g).toString(16);
  if ( g.length == 1 ) {
    g = "0" + g;
  }

  b = Math.floor(b).toString(16);
  if ( b.length == 1 ) {
    b = "0" + b;
  }

  return "#" + r + g + b;
}



function randhuecolor(cmin, cmax)
{
  var x = Math.random() * 6;
  var i = Math.floor( x ), r = 0, g = 0, b = 0;

  if ( cmin === undefined || cmin === null ) {
    cmin = 0;
  }
  if ( cmax === undefined || cmax === null ) {
    cmax = 255;
  }
  var cvar = cmax - cmin + 1;
  x -= i;
  if ( i < 1 ) { // red to yellow
    r = cmax;
    g = cmin + Math.floor( cvar * x );
  } else if ( i < 2 ) { // yellow to green
    r = cmin + Math.floor( cvar * (1 - x) );
    g = cmax;
  } else if ( i < 3 ) { // green to cyan
    g = cmax;
    b = cmin + Math.floor( cvar * x );
  } else if ( i < 4 ) { // cyan to blue
    g = cmin + Math.floor( cvar * (1 - x) );
    b = cmax;
  } else if ( i < 5 ) { // blue to magenta
    b = cmax;
    r = cmin + Math.floor( cvar * x );
  } else {
    b = cmin + Math.floor( cvar * (1 - x) );
    r = cmax;
  }
  return rgb2str(r, g, b);
}



/* darken a color */
function darkenColor(colorStr) {
  // Defined in dygraph-utils.js
  var color = Dygraph.toRGB_(colorStr);
  color.r = Math.floor((255 + color.r) / 2);
  color.g = Math.floor((255 + color.g) / 2);
  color.b = Math.floor((255 + color.b) / 2);
  return 'rgb(' + color.r + ',' + color.g + ',' + color.b + ')';
}



/* bar chart for Dygraph */
function barChartPlotter(e) {
  var ctx = e.drawingContext;
  var points = e.points;
  var y_bottom = e.dygraph.toDomYCoord(0);

  ctx.fillStyle = darkenColor(e.color);

  // Find the minimum separation between x-values.
  // This determines the bar width.
  var min_sep = Infinity;
  for (var i = 1; i < points.length; i++) {
    var sep = points[i].canvasx - points[i - 1].canvasx;
    if (sep < min_sep) min_sep = sep;
  }
  var bar_width = Math.floor(2.0 / 3 * min_sep);

  // Do the actual plotting.
  for (var i = 0; i < points.length; i++) {
    var p = points[i];
    var center_x = p.canvasx;

    ctx.fillRect(center_x - bar_width / 2, p.canvasy,
        bar_width, y_bottom - p.canvasy);

    ctx.strokeRect(center_x - bar_width / 2, p.canvasy,
        bar_width, y_bottom - p.canvasy);
  }
}
