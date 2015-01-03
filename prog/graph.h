#ifndef GRAPH_H__
#define GRAPH_H__



typedef struct {
  int n;
  char *mat;
  int *cid;
  int *queue;
  int *chist; /* histogram by the cluster size */
} graph_t;



static graph_t *graph_open(int n)
{
  graph_t *g;

  xnew(g, 1);
  g->n = n;
  xnew(g->mat, n * n);
  xnew(g->cid, n);
  xnew(g->queue, n);
  xnew(g->chist, n + 1);
  return g;
}



static void graph_close(graph_t *g)
{
  free(g->mat);
  free(g->cid);
  free(g->queue);
  free(g->chist);
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
    g->chist[end] += 1;
  }
}



/* clear the cluster histogram */
__inline static void graph_chist_clear(graph_t *g)
{
  int i;

  for ( i = 0; i < g->n; i++ ) g->chist[i] = 0;
}



/* print the cluster histogram */
__inline static void graph_chist_print(const graph_t *g)
{
  int i;

  for ( i = 1; i < g->n; i++ )
    if ( g->chist[i] )
      printf("%4d: %d\n", i, g->chist[i]);
}



#endif /* GRAPH_H__ */

