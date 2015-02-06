#ifndef LJCLS_H__
#define LJCLS_H__



/* clustering based energy */



#include "ljcore.h"
#include "graph.h"
#include "hmc.h"
#include "wl.h"



typedef struct {
  lj_t *lj;

  graph_t *g, *g2;
  double rcls; /* cluster cutoff */
  const double *vcls; /* bias potential */
  int cseed; /* seed of cluster */
  double (*x2)[D];
} ljcls_t;




__inline static ljcls_t *ljcls_open(lj_t *lj,
    double rcls, const double *vcls)
{
  ljcls_t *c;
  int n = lj->n;

  xnew(c, 1);
  c->lj = lj;

  c->g = graph_open(n);
  c->g2 = graph_open(n);
  c->rcls = rcls;
  c->vcls = vcls;
  c->cseed = 0;
  xnew(c->x2, n);
  return c;
}



__inline static void ljcls_close(ljcls_t *c)
{
  graph_close(c->g);
  graph_close(c->g2);
  free(c->x2);
  free(c);
}



#define ljcls_mkgraph(c, g) ljcls_mkgraph_low(c->lj, g, c->rcls)

/* build a graph and do clustering */
static void ljcls_mkgraph_low(lj_t *lj, graph_t *g, double rm)
{
  int i, j, n = lj->n;
  double rm2 = rm * rm;

  graph_empty(g);
  for ( i = 0; i < n; i++ )
    for ( j = i + 1; j < n; j++ )
      if ( lj->r2ij[i*n + j] < rm2 )
        graph_link(g, i, j);
  graph_clus(g);
}



#define ljcls_mkgraph2(c, g, k) ljcls_mkgraph2_low(c->lj, g, k, c->rcls)

/* build a graph with the distances from k computed r2i */
static void ljcls_mkgraph2_low(lj_t *lj, graph_t *g, int k, double rm)
{
  int i, j, n = lj->n;
  double r2, rm2 = rm * rm;

  graph_empty(g);
  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      if ( i == k ) {
        r2 = lj->r2i[j];
      } else if ( j == k ) {
        r2 = lj->r2i[i];
      } else {
        r2 = lj->r2ij[i*n + j];
      }
      if ( r2 < rm2 )
        graph_link(g, i, j);
    }
  }
  graph_clus(g);
}



/* compute the cluster energy
 * call ljcls_mkgraph() first */
__inline static double ljcls_eclus(const ljcls_t *c, const graph_t *g)
{
  return c->vcls ? c->vcls[ graph_getcsize(g, c->cseed) ] : 0.0;
}



/* Metropolis algorithm */
__inline static int ljcls_metro(ljcls_t *c, double amp, double bet)
{
  lj_t *lj = c->lj;
  int i, acc = 0;
  double xi[D], r, du = 0, dutot, dvir = 0;
  double ucls, ducls = 0;

  i = lj_randmv(lj, xi, amp);
  du = lj_depot(lj, i, xi, &dvir);

  ljcls_mkgraph2(c, c->g2, i);

  /* compute the clustering energy */
  ucls = ljcls_eclus(c, c->g2);
  /* since the clustering potential might have been changed
   * we have to recompute the old clustering energy */
  ducls = ucls - ljcls_eclus(c, c->g);

  dutot = bet * du + ducls;
  if ( dutot < 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp( -dutot ) );
  }
  if ( acc ) {
    lj_commit(lj, i, xi, du, dvir);
    graph_copy(c->g, c->g2);
    return 1;
  }
  return 0;
}



/* attempt to change the seed for clustering */
__inline static int ljcls_changeseed(ljcls_t *c, const graph_t *g)
{
  lj_t *lj = c->lj;
  int i, sz0, sz1, n = lj->n, acc = 0;

  i = (c->cseed + 1 + (int) (rand01() * (n - 1))) % n;
  sz0 = graph_getcsize(g, c->cseed);
  sz1 = graph_getcsize(g, i);
  if ( sz0 == sz1 || c->vcls == NULL ) {
    acc = 1;
  } else {
    double dv = c->vcls[ sz1 ] - c->vcls[ sz0 ];
    if ( dv < 0 ) {
      acc = 1;
    } else {
      double r = rand01();
      acc = (r < exp( -dv ));
    }
  }
  if ( acc ) {
    c->cseed = i;
  }
  return acc;
}



/* a step of HMC
 * `*csize1` gives the current cluster size on return */
__inline static int ljcls_hmc(ljcls_t *c, hmc_t *hmc, int *csize1)
{
  lj_t *lj = c->lj;
  /* compute the current cluster size */
  int csize = graph_getcsize(c->g, c->cseed);
  int acc, idat[2] = { csize, c->cseed };
  /* hmc->idat[0] is the previous size */
  double dv = 0;

  if ( c->vcls != NULL ) {
    dv = c->vcls[ csize ] - c->vcls[ hmc->idat[0] ];
  }

  if ( dv <= 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = ( r < exp( -dv ) );
  }
  if ( acc ) {
    hmc_push(hmc, lj->x, lj->v, lj->f, idat, &lj->epot);
  } else {
    hmc_pop(hmc, lj->x, lj->v, lj->f, idat, &lj->epot, 1);
    c->cseed = idat[1];
    /* lj->r2ij should be refreshed in the next step */
    //lj_force(lj);
    //lj_mkgraph(lj, c->g);
  }
  *csize1 = idat[0];
  return acc;
}



/* wrap coordinates such that particles
 * in the same cluster stay close */
__inline static int ljcls_wrapclus(ljcls_t *c,
    double (*xin)[D], double (*xout)[D], graph_t *g)
{
  int ic, i, j, n = c->lj->n, head, end;
  double l = c->lj->l, invl = 1 / l, dx[D];

  ljcls_mkgraph(c, g);
  for ( i = 0; i < n; i++ ) {
    vwrap( vcopy(xout[i], xin[i]), l );
  }

  for ( ic = 0; ic < g->nc; ic++ ) {
    /* find the seed of this cluster */
    for ( i = 0; i < n; i++ )
      if ( g->cid[i] == ic )
        break;
    if ( i >= n ) {
      fprintf(stderr, "no particle belongs to cluster %d\n", ic);
      return -1;
    }

    g->queue[ head = 0 ] = i;
    g->cid[ i ] = -1;
    end = 1;
    for ( ; head < end; head++ ) {
      i = g->queue[head];
      /* add neighbors of i into the queue */
      for ( j = 0; j < n; j++ ) {
        if ( graph_linked(g, i, j) && g->cid[j] >= 0 ) {
          g->cid[j] = -1;
          g->queue[ end++ ] = j;
          /* align j with i */
          lj_vpbc(vdiff(dx, xout[j], xout[i]), l, invl);
          vinc( vcopy(xout[j], xout[i]), dx );
        }
      }
    }
    if ( end != g->csize[ic] ) {
      fprintf(stderr, "cluster %d: size %d vs %d\n", ic, end, g->csize[ic]);
    }
  }
  return 0;
}



__inline static int ljcls_writepos(ljcls_t *c,
    double (*x)[D], double (*v)[D],
    const char *fn, int wrap)
{
  lj_t *lj = c->lj;

  if ( wrap ) {
    ljcls_wrapclus(c, x, c->x2, c->g2);
  } else {
    lj_wrapbox(lj, x, c->x2);
  }

  return lj_writepos(lj, c->x2, v, fn);
}



/* randomly swap the velocities of k pairs of particles */
#define lj_vscramble(lj, v, k) md_vscramble(v, NULL, lj->n, k)



#endif /* LJCLS_H__ */

