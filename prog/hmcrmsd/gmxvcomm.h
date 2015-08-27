#ifndef GMXVCOMM_H__
#define GMXVCOMM_H__


#include "gmx_ga2la.h"

/* the following macros are defined in domdec.c */
#ifndef DDRANK
#define DDRANK(dd, rank)    (rank)
#endif
#ifndef DDMASTERRANK
#define DDMASTERRANK(dd)   (dd->masterrank)
#endif



/* communicate vectors of a set of special atoms */
typedef struct {
  t_commrec *cr;
  t_state *state_local, *state_global;
  int n;
  int *index; /* global indices of the special atoms */
  int lcnt; /* number of special atoms in this node */
  int *lcnt_m; /* master collector of `lcnt` */
  int *lwho; /* lwho[i] is gives the index of special atoms
                `0 <= i < lcnt` */
  int *lwho_m; /* master collector of `lwho` */
  int *la; /* la[i] is gives the local index of the ith
              special atom in the buffer, `0 <= i < lcnt` */
  int *ga; /* ga[i] is gives the global index of the ith
              special atom in the buffer, `0 <= i < lcnt` */
  rvec *lx; /* local vectors */
  rvec *lx_m; /* master collector of `lx` */
} gmxvcomm_t;



static gmxvcomm_t *gmxvcomm_open(int n, int *arr,
    t_commrec *cr, t_state *state_local, t_state *state_global)
{
  gmxvcomm_t *g;
  int i, sz = 0;

  snew(g, 1);
  g->n = n;
  snew(g->index, n);
  for ( i = 0; i < n; i++ ) {
    g->index[i] = arr[i];
  }
  g->cr = cr;
  g->state_local = state_local;
  g->state_global = state_global;
  snew(g->lwho, n);
  snew(g->la, n);
  snew(g->ga, n);
  snew(g->lx, n);

  if ( MASTER(cr) ) {
    if ( DOMAINDECOMP(cr) ) {
      MPI_Comm_size(cr->mpi_comm_mygroup, &sz);
    } else {
      sz = 1;
    }
    snew(g->lcnt_m, sz);
    /* n should be enough */
    snew(g->lwho_m, n);
    snew(g->lx_m, n);
  }
  return g;
}



static void gmxvcomm_close(gmxvcomm_t *g)
{
  free(g->index);
  free(g->lcnt_m);
  free(g->lwho);
  free(g->la);
  free(g->ga);
  free(g->lx);
  free(g);
}



static int gmxvcomm_print(gmxvcomm_t *g)
{
  int i, ig, il;
  t_commrec *cr = g->cr;
  rvec *xl = g->state_local->x;
  rvec *xg = g->state_global->x;

  if ( !DOMAINDECOMP(cr) ) {
    if ( MASTER(cr) ) {
      for ( i = 0; i < g->n; i++ ) {
        ig = g->index[i];
        fprintf(stderr, "index %d, %g %g %g\n", ig, xg[ig][0], xg[ig][1], xg[ig][2]);
      }
    }
    return 0;
  }

  for ( i = 0; i < g->n; i++ ) {
    ig = g->index[i];
    if ( ga2la_get_home(cr->dd->ga2la, ig, &il) ) {
      fprintf(stderr, "%d/%d: index %d, lwho %d, %g %g %g\n",
          cr->nodeid, cr->nnodes, ig, il,
          xl[il][0], xl[il][1], xl[il][2]);
    }
    if ( DDMASTER(cr->dd) ) {
      fprintf(stderr, "index %d, %g %g %g\n", ig, xg[ig][0], xg[ig][1], xg[ig][2]);
    }
  }

  return 0;
}



/* let the master know how many special atoms are in each node
 * cf. dd_collect_vec_sendrecv() defined in domdec.c */
static int gmxvcomm_gatherid(gmxvcomm_t *g)
{
  int i, ig, il, iw, lcnt;
  t_commrec *cr = g->cr;
  gmx_domdec_t *dd = cr->dd;

  if ( !DOMAINDECOMP(cr) ) return 0;

  g->lcnt = 0;
  for ( i = 0; i < g->n; i++ ) {
    ig = g->index[i];
    if ( ga2la_get_home(cr->dd->ga2la, ig, &il) ) {
      g->lwho[ g->lcnt ] = i;
      g->la[ g->lcnt ] = il;
      g->ga[ g->lcnt ] = ig;
      g->lcnt++;
    }
  }

#ifdef GMX_MPI
  MPI_Gather(&(g->lcnt), 1, MPI_INT, g->lcnt_m, 1, MPI_INT,
      DDMASTERRANK(dd), dd->mpi_comm_all);
#endif

  if ( DDMASTER(dd) ) {
    /* master collects `lwho` to `lwho_m` */
    for ( iw = 0, i = 0; i < dd->nnodes; i++ ) {
      if ( i != dd->rank && (lcnt = g->lcnt_m[i]) > 0 ) {
#ifdef GMX_MPI
        MPI_Recv(g->lwho_m + iw, lcnt, MPI_INT,
            DDRANK(dd, i), i, dd->mpi_comm_all, MPI_STATUS_IGNORE);
#endif
        //printf("i %d, cnt %d, iw %d\n", i, lcnt, iw);
        //if ( lcnt > 0 ) { int j; for ( j = 0; j < lcnt; j++ ) printf(" %d", g->lwho[iw + j]); printf("\n"); }
        iw += lcnt;
      }
    }
  } else {
    /* slaves send `lwho` */
    if ( g->lcnt > 0 ) {
#ifdef GMX_MPI
      MPI_Send(g->lwho, g->lcnt, MPI_INT,
          DDMASTERRANK(dd), dd->rank, dd->mpi_comm_all);
#endif
    }
  }

  return 0;
}



/* collect vectors (x),
 * must call `gmxvcomm_gatherid()` first */
static int gmxvcomm_gatherv(gmxvcomm_t *g)
{
  int i, iw, lcnt;
  t_commrec *cr = g->cr;
  gmx_domdec_t *dd = cr->dd;

  if ( DDMASTER(dd) ) {
    /* master collects the vectors */
    for ( iw = 0, i = 0; i < dd->nnodes; i++ ) {
      if ( i != dd->rank && (lcnt = g->lcnt_m[i]) > 0 ) {
#ifdef MPI
        MPI_Recv(g->lx_m + iw, lcnt * sizeof(rvec), MPI_BYTE,
            DDRANK(dd, i), i, dd->mpi_comm_all, MPI_STATUS_IGNORE);
#endif
        iw += lcnt;
      }
    }
  } else {
    /* slaves send `lx` */
    if ( g->lcnt > 0 ) {
#ifdef MPI
      MPI_Send(g->lx, g->lcnt * sizeof(rvec), MPI_BYTE,
          DDMASTERRANK(dd), dd->rank, dd->mpi_comm_all);
#endif
    }
  }

  return 0;
}



/* distribute vectors (f) */
static int gmxvcomm_scatterv(gmxvcomm_t *g)
{
  int i, iw, lcnt;
  t_commrec *cr = g->cr;
  gmx_domdec_t *dd = cr->dd;

  if ( DDMASTER(dd) ) {
    /* master collects the vectors */
    for ( iw = 0, i = 0; i < dd->nnodes; i++ ) {
      if ( i != dd->rank && (lcnt = g->lcnt_m[i]) > 0 ) {
#ifdef MPI
        MPI_Send(g->lx_m + iw, lcnt * sizeof(rvec), MPI_BYTE,
            DDRANK(dd, i), i, dd->mpi_comm_all);
#endif
        iw += lcnt;
      }
    }
  } else {
    /* slaves send `lx` */
    if ( g->lcnt > 0 ) {
#ifdef MPI
      MPI_Recv(g->lx, g->lcnt * sizeof(rvec), MPI_BYTE,
          DDMASTERRANK(dd), dd->rank, dd->mpi_comm_all, MPI_STATUS_IGNORE);
#endif
    }
  }

  return 0;
}



#endif /* defined(GMXVCOMM__) */
