#ifndef GRAPH_H__
#define GRAPH_H__



typedef struct {
  int n;
  char *mat;
  int nc; /* number of connected compoments */
  int *csize; /* csize[i] is the number of particles in component i */
  int *cid; /* cid[i] is the component id of particle i */
  int *queue; /* queue used in finding the connected components */
} graph_t;



static graph_t *graph_open(int n)
{
  graph_t *g;

  xnew(g, 1);
  g->n = n;
  xnew(g->mat, n * n);
  xnew(g->csize, n);
  xnew(g->cid, n);
  xnew(g->queue, n);
  return g;
}



static void graph_close(graph_t *g)
{
  free(g->mat);
  free(g->csize);
  free(g->cid);
  free(g->queue);
  free(g);
}



static void graph_empty(graph_t *g)
{
  memset(g->mat, 0, sizeof(g->mat[0]) * g->n * g->n);
}



static void graph_link(graph_t *g, int i, int j)
{
  g->mat[i*g->n + j] = '\1';
  g->mat[j*g->n + i] = '\1';
}



static int graph_linked(const graph_t *g, int i, int j)
{
  return g->mat[i*g->n + j];
}


__inline static int graph_copy(graph_t *g, const graph_t *g2)
{
  int i, n = g2->n;

  if ( g->n != n ) {
    fprintf(stderr, "cannot copy graphs of different sizes %d vs %d\n", g->n, n);
    return -1;
  }
  for ( i = 0; i < n * n; i++ )
    g->mat[i] = g2->mat[i];

  /* copy the cluster information */
  g->nc = g2->nc;
  for ( i = 0; i < n; i++ ) {
    g->csize[i] = g2->csize[i];
    g->cid[i] = g2->cid[i];
  }
  return 0;
}


/* do clustering */
__inline static void graph_clus(graph_t *g)
{
  int i, j, n = g->n, ic;
  int head, end;

  for ( i = 0; i < n; i++ ) g->cid[i] = -1;

  for ( ic = 0; ; ic++ ) {
    /* find the first free atom */
    for ( i = 0; i < n; i++ )
      if ( g->cid[i] < 0 )
        break;
    if ( i >= n ) break;

    /* grow a cluster starting from i */
    g->queue[ head = 0 ] = i;
    g->cid[i] = ic;
    end = 1;
    for ( ; head < end; head++ ) {
      i = g->queue[head];
      /* add neighbors of i in to the queue */
      for ( j = 0; j < n; j++ ) {
        if ( graph_linked(g, i, j) && g->cid[j] < 0 ) {
          g->cid[j] = ic;
          g->queue[ end++ ] = j;
        }
      }
    }
    /* now `end` is the number of particles in cluster `ic` */
    g->csize[ic] = end;
  }
  g->nc = ic;
}



__inline static void graph_clus_print(const graph_t *g)
{
  int i;
  printf("%d clusters: ", g->nc);
  for ( i = 0; i < g->nc; i++ )
    printf("%d ", g->csize[i]);
  printf("\n");
}



#endif /* GRAPH_H__ */

