#ifndef LJMIXCLS_H__
#define LJMIXCLS_H__



/* clustering based energy */



#include "ljmixcore.h"
#include "graph.h"
#include "hmc.h"
#include "wl.h"



typedef struct {
  ljmix_t *lj;

  graph_t *g, *g2;
  double rcls; /* cluster cutoff */
  const double *vcls; /* bias potential */
  int cseed; /* seed of cluster */
  double (*x2)[D];
} ljmixcls_t;




__inline static ljmixcls_t *ljmixcls_open(ljmix_t *lj,
    double rcls, const double *vcls)
{
  ljmixcls_t *c;
  int np = lj->np[0];

  xnew(c, 1);
  c->lj = lj;

  c->g = graph_open(np);
  c->g2 = graph_open(np);
  c->rcls = rcls;
  c->vcls = vcls;
  c->cseed = 0;
  xnew(c->x2, lj->n);
  return c;
}



__inline static void ljmixcls_close(ljmixcls_t *c)
{
  graph_close(c->g);
  graph_close(c->g2);
  free(c->x2);
  free(c);
}



#define ljmixcls_mkgraph(c, g) ljmixcls_mkgraph_low(c->lj, g, c->rcls)

/* build a graph and do clustering */
static void ljmixcls_mkgraph_low(ljmix_t *lj, graph_t *g, double rm)
{
  int i, j, np = lj->np[0];
  double rm2 = rm * rm;

  graph_empty(g);
  for ( i = 0; i < np; i++ ) {
    for ( j = i + 1; j < np; j++ ) {
      if ( lj->r2ij[i*np + j] < rm2 ) {
        graph_link(g, i, j);
      }
    }
  }
  graph_clus(g);
}



#define ljmixcls_mkgraph2(c, g, k) ljmixcls_mkgraph2_low(c->lj, g, k, c->rcls)

/* build a graph with the distances from k computed r2i */
static void ljmixcls_mkgraph2_low(ljmix_t *lj, graph_t *g, int k, double rm)
{
  int i, j, np = lj->np[0];
  double r2, rm2 = rm * rm;

  graph_empty(g);
  for ( i = 0; i < np; i++ ) {
    for ( j = i + 1; j < np; j++ ) {
      if ( i == k ) {
        r2 = lj->r2i[j];
      } else if ( j == k ) {
        r2 = lj->r2i[i];
      } else {
        r2 = lj->r2ij[i*np + j];
      }
      if ( r2 < rm2 ) {
        graph_link(g, i, j);
      }
    }
  }
  graph_clus(g);
}



/* compute the cluster energy
 * call ljmixcls_mkgraph() first */
__inline static double ljmixcls_eclus(const ljmixcls_t *c, const graph_t *g)
{
  return c->vcls ? c->vcls[ graph_getcsize(g, c->cseed) ] : 0.0;
}



/* Metropolis algorithm */
__inline static int ljmixcls_metro(ljmixcls_t *c, double amp, double bet)
{
  ljmix_t *lj = c->lj;
  int i, acc = 0;
  double xi[D], r, du, dutot, dvir = 0;
  double ucls, ducls = 0;

  i = ljmix_randmv(lj, xi, amp);
  du = ljmix_depot(lj, i, xi, &dvir);

  ljmixcls_mkgraph2(c, c->g2, i);

  /* compute the clustering energy */
  ucls = ljmixcls_eclus(c, c->g2);
  /* since the clustering potential might have been changed
   * we have to recompute the old clustering energy */
  ducls = ucls - ljmixcls_eclus(c, c->g);

  dutot = bet * du + ducls;
  if ( dutot < 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp( -dutot ) );
  }
  if ( acc ) {
    ljmix_commit(lj, i, xi, du, dvir);
    graph_copy(c->g, c->g2);
    return 1;
  }
  return 0;
}



/* attempt to change the seed for clustering */
__inline static int ljmixcls_changeseed(ljmixcls_t *c, const graph_t *g)
{
  ljmix_t *lj = c->lj;
  int i, sz0, sz1, np = lj->np[0], acc = 0;

  i = (c->cseed + 1 + (int) (rand01() * (np - 1))) % np;
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
__inline static int ljmixcls_hmc(ljmixcls_t *c, hmc_t *hmc, int *csize1)
{
  ljmix_t *lj = c->lj;
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
__inline static int ljmixcls_wrapclus(ljmixcls_t *c,
    double (*xin)[D], double (*xout)[D], graph_t *g)
{
  int ic, i, j, n = c->lj->n, np = c->lj->np[0], head, end;
  double l = c->lj->l, invl = 1 / l, dx[D];

  ljmixcls_mkgraph(c, g);
  for ( i = 0; i < n; i++ ) {
    vwrap( vcopy(xout[i], xin[i]), l );
  }

  for ( ic = 0; ic < g->nc; ic++ ) {
    /* find the seed of this cluster */
    for ( i = 0; i < np; i++ )
      if ( g->cid[i] == ic )
        break;
    if ( i >= np ) {
      fprintf(stderr, "no particle belongs to cluster %d\n", ic);
      return -1;
    }

    g->queue[ head = 0 ] = i;
    g->cid[ i ] = -1;
    end = 1;
    for ( ; head < end; head++ ) {
      i = g->queue[head];
      /* add neighbors of i into the queue */
      for ( j = 0; j < np; j++ ) {
        if ( graph_linked(g, i, j) && g->cid[j] >= 0 ) {
          g->cid[j] = -1;
          g->queue[ end++ ] = j;
          /* align j with i */
          ljmix_vpbc(vdiff(dx, xout[j], xout[i]), l, invl);
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



__inline static int ljmixcls_writepos(ljmixcls_t *c,
    double (*x)[D], double (*v)[D],
    const char *fn, int wrap)
{
  ljmix_t *lj = c->lj;

  if ( wrap ) {
    ljmixcls_wrapclus(c, x, c->x2, c->g2);
  } else {
    ljmix_wrapbox(lj, x, c->x2);
  }

  return ljmix_writepos(lj, c->x2, v, fn);
}



/* randomly swap the velocities of k pairs of particles */
#define ljmix_vscramble(lj, v, k) md_vscramble(v, NULL, lj->n, k)



#endif /* LJMIXCLS_H__ */

