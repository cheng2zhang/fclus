#ifndef GMXVCOMM_H__
#define GMXVCOMM_H__



/* Communicator for a set of coordinates */



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
  double (*lx)[3]; /* local vectors */
  double (*lx_m)[3]; /* master collector of `lx` */
} gmxvcomm_t;



static gmxvcomm_t *gmxvcomm_open(int n, int *arr,
    t_commrec *cr)
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
  snew(g->lwho, n);
  snew(g->la, n);
  snew(g->ga, n);
  snew(g->lx, n);

  if ( MASTER(cr) ) {
    if ( DOMAINDECOMP(cr) ) {
#ifdef GMX_MPI
      MPI_Comm_size(cr->mpi_comm_mygroup, &sz);
#endif
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
  sfree(g->index);
  sfree(g->lcnt_m);
  sfree(g->lwho);
  sfree(g->la);
  sfree(g->ga);
  sfree(g->lx);
  sfree(g);
}



/* let the master know how many special atoms are in each node
 * cf. dd_collect_vec_sendrecv() defined in domdec.c */
static int gmxvcomm_gatherid(gmxvcomm_t *g)
{
  int i, ig, il, iw, lcnt;
  t_commrec *cr = g->cr;
  gmx_domdec_t *dd = cr->dd;

  if ( !DOMAINDECOMP(cr) ) return 0;

  /* locally compute `lwho`, which is a list
   * collecting special atoms on the node */
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
      if ( (lcnt = g->lcnt_m[i]) > 0 ) {
        if ( i != dd->rank ) {
#ifdef GMX_MPI
          MPI_Recv(g->lwho_m + iw, lcnt, MPI_INT,
              DDRANK(dd, i), i, dd->mpi_comm_all, MPI_STATUS_IGNORE);
#endif
        } else {
          memcpy(g->lwho_m + iw, g->lwho, lcnt * sizeof(int));
        }
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
  gmx_domdec_t *dd;

  if ( !DOMAINDECOMP(cr) ) return 0;

  dd = cr->dd;
  if ( DDMASTER(dd) ) {
    /* master collects the vectors */
    for ( iw = 0, i = 0; i < dd->nnodes; i++ ) {
      if ( (lcnt = g->lcnt_m[i]) > 0 ) {
        if ( i != dd->rank ) {
#ifdef GMX_MPI
          MPI_Recv(g->lx_m + iw, lcnt * 3, MPI_DOUBLE,
              DDRANK(dd, i), i, dd->mpi_comm_all, MPI_STATUS_IGNORE);
#endif
        } else {
          memcpy(g->lx_m + iw, g->lx, lcnt * sizeof(g->lx[0]));
        }
        iw += lcnt;
      }
    }
  } else {
    /* slaves send `lx` */
    if ( g->lcnt > 0 ) {
#ifdef GMX_MPI
      MPI_Send(g->lx, g->lcnt * 3, MPI_DOUBLE,
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
  gmx_domdec_t *dd;

  if ( !DOMAINDECOMP(cr) ) return 0;

  dd = cr->dd;
  if ( DDMASTER(dd) ) {
    /* master collects the vectors */
    for ( iw = 0, i = 0; i < dd->nnodes; i++ ) {
      if ( (lcnt = g->lcnt_m[i]) > 0 ) {
        if ( i != dd->rank ) {
#ifdef GMX_MPI
          MPI_Send(g->lx_m + iw, lcnt * 3, MPI_DOUBLE,
              DDRANK(dd, i), i, dd->mpi_comm_all);
#endif
        } else {
          memcpy(g->lx, g->lx_m + iw, lcnt * sizeof(g->lx[0]));
        }
        iw += lcnt;
      }
    }
  } else {
    /* slaves send `lx` */
    if ( g->lcnt > 0 ) {
#ifdef GMX_MPI
      MPI_Recv(g->lx, g->lcnt * 3, MPI_DOUBLE,
          DDMASTERRANK(dd), dd->rank, dd->mpi_comm_all, MPI_STATUS_IGNORE);
#endif
    }
  }

  return 0;
}



/* find the rank of special atom i
 * only call this on the master */
static int gmxvcomm_where(gmxvcomm_t *g, int id)
{
  int i, j, iw, lcnt;
  t_commrec *cr = g->cr;
  gmx_domdec_t *dd;

  if ( !DOMAINDECOMP(cr) ) return 0;

  dd = cr->dd;
  for ( iw = 0, i = 0; i < dd->nnodes; i++ ) {
    if ( (lcnt = g->lcnt_m[i]) > 0 ) {
      for ( j = 0; j < lcnt; j++ ) {
        if ( g->lwho_m[iw + j] == id )
          return i;
      }
      iw += lcnt;
    }
  }

  return 0;
}



#endif /* defined(GMXVCOMM_H__) */

