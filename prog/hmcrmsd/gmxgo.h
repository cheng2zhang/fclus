#ifndef GMXGO_H__
#define GMXGO_H__



/* Go model */
#define D 3
#include "mat.h"
#include "wl.h"
#include "mtrand.h"
#include "gmxvcomm.h"
#include "gmxgomodel.h"



typedef struct {
  int n; /* number of atoms */
  int *index; /* index[0..n-1]: global indices of the Go atoms */
  double *mass; /* mass[0..n-1] mass */
  double masstot;
  double (*x)[3]; /* x[0..n-1] */
  double (*f)[3]; /* f[0..n-1] bias force */
  double (*xref)[3];
  double (*xf)[3]; /* fit structure */
  double (*x1)[3]; /* x1[0..n-1] */
  double (*x2)[3]; /* x2[0..n-1] */
  gmxvcomm_t *gvc;
  gmxgomodel_t model[1];
  double kT;
  wl_t *wl;
  double *rhis; /* histogram */
  double rhis_dx, rhis_max; /* spacing and maximal */
  int rhis_n;
  t_commrec *cr;
  double dvdx;

  /* HMC state variables */
  int stn, stncap;
  rvec *stx;
  rvec *stv;
  rvec *stf;
  /*
  matrix box;
  */

  /* string buffer to print step */
  char sbuf[STEPSTRSIZE];
} gmxgo_t;



/* extract the coordinates of the special atoms */
static int gmxgo_build_master(gmxgo_t *go, gmx_mtop_t *mtop)
{
  int itp, ib, id, ia, ncap, i;
  gmx_moltype_t *mt;
  gmx_molblock_t *mb = mtop->molblock;

  /* search for the moltype named "Protein" */
  for ( itp = 0; itp < mtop->nmoltype; itp++ ) {
    if ( strncmp(mtop->moltype[itp].name[0], "Protein", 7) == 0 ) {
      break;
    }
  }
  /* assume the default moltype to 0 */
  if ( itp >= mtop->nmoltype ) {
    itp = 0;
  }
  mt = mtop->moltype + itp;

  /* find the molblock correpsonding to the moltype */
  id = 0;
  for ( ib = 0; ib < mtop->nmolblock; ib++ ) {
    mb = mtop->molblock + ib;
    if ( mb->type == itp ) break;
    id += mb->nmol * mb->natoms_mol;
  }

#define BLKSZ 8
  /* search CA atoms from the moltype */
  go->n = 0;
  go->masstot = 0;
  ncap = BLKSZ;
  snew(go->index, ncap);
  snew(go->mass, ncap);
  for ( ia = 0; ia < mt->atoms.nr; ia++ ) {
    char *atnm = mt->atoms.atomname[ia][0];
    int sel = 0;

    if ( strcmp(atnm, "CA") == 0 ) {
      sel = 1;
    }

    if ( sel ) {
      /* reallocate the memory */
      if ( go->n >= ncap ) {
        ncap += BLKSZ;
        srenew(go->index, ncap);
        srenew(go->mass, ncap);
      }

      go->index[go->n] = id + ia;
      go->mass[go->n] = mt->atoms.atom[ia].m;
      go->masstot += go->mass[go->n];
      go->n += 1;
    }
  }

  fprintf(stderr, "%d spatoms:", go->n);
  for ( i = 0; i < go->n; i++ ) {
    fprintf(stderr, " %d(%g)", go->index[i], go->mass[i]);
  }
  fprintf(stderr, "; masstot %g\n", go->masstot);

  return 0;
}



/* extract the coordinates of the special atoms */
static int gmxgo_build(gmxgo_t *go, gmx_mtop_t *mtop, t_commrec *cr)
{
  if ( MASTER(cr) ) {
    gmxgo_build_master(go, mtop);
  }

  if ( PAR(cr) ) {
    /* Note: gmx_bcast() works for cr->mpi_comm_mygroup
     * see gmxlib/network.c */
    gmx_bcast(sizeof(int), &go->n, cr);
    if ( !MASTER(cr) ) {
      snew(go->index, go->n);
      snew(go->mass, go->n);
    }
    gmx_bcast(sizeof(go->index[0]) * go->n, go->index, cr);
    gmx_bcast(sizeof(go->mass[0]) * go->n, go->mass, cr);
    gmx_bcast(sizeof(go->masstot), &go->masstot, cr);
  }

  return 0;
}



/* load the reference coordinates */
static int gmxvcomm_loadxref(gmxgo_t *go, const char *fnpdb)
{
  FILE *fp;
  char s[256];
  int id = 0;
  float x[3];

  if ( (fp = fopen(fnpdb, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fnpdb);
    return -1;
  }

  while ( fgets(s, sizeof s, fp) ) {
    /* a model is finished, we implicitly include ENDMDL */
    if ( strncmp(s, "TER", 3) == 0 || strncmp(s, "END", 3) == 0 )
      break;

    if ( strncmp(s, "ATOM ", 5) != 0 )
      continue;

    if ( s[16] != ' ' && s[16] != 'A' )
      continue;

    if ( strncmp(s + 12, " CA ", 4) != 0 )
      continue;

    if ( id >= go->n ) {
      fprintf(stderr, "Warning: additional atoms exists in %s\n", fnpdb);
      break;
    }

    /* scan the coordinates */
    sscanf(s + 30 , "%f%f%f", x, x + 1, x + 2);

    /* convert the unit from angstroms to nanometers */
    go->xref[id][0] = x[0] / 10;
    go->xref[id][1] = x[1] / 10;
    go->xref[id][2] = x[2] / 10;
    id += 1;
  }

  if ( id < go->n ) {
    fprintf(stderr, "Warning: missing atoms %s, %d < %d\n", fnpdb, id, go->n);
    return -1;
  }

  return 0;
}



/* initialize */
static void gmxgo_inithmc(gmxgo_t *go, gmx_mtop_t *mtop)
{
  int ib, n = 0;
  gmx_molblock_t *mb = mtop->molblock;
  t_commrec *cr = go->cr;

  /* count the number of atoms */
  for ( ib = 0; ib < mtop->nmolblock; ib++ ) {
    n += mb[ib].nmol * mb[ib].natoms_mol;
  }

  if ( DOMAINDECOMP(cr) ) {
    /* we allow some margin */
    n = (int) ((n + 3 * sqrt(n)) / cr->nnodes);
  }

  go->stncap = n;
  go->stn = 0;
  snew(go->stx, n);
  snew(go->stv, n);
  snew(go->stf, n);
}



/* tp: the temperature */
static gmxgo_t *gmxgo_open(gmx_mtop_t *mtop, t_commrec *cr,
    double tp, const char *fncfg)
{
  gmxgo_t *go;

  snew(go, 1);

  go->cr = cr;
  gmxgo_build(go, mtop, cr);
  if ( PAR(cr) ) {
    go->gvc = gmxvcomm_open(go->n, go->index, cr);
  } else {
    go->gvc = NULL;
  }

  snew(go->x, go->n);
  snew(go->f, go->n);
  snew(go->xref, go->n);
  snew(go->xf, go->n);
  snew(go->x1, go->n);
  snew(go->x2, go->n);

  /* load the reference */
  if ( MASTER(cr) ) {
    gmxgomodel_t *m = go->model;

    gmxgomodel_default(m);
    gmxgomodel_load(m, fncfg);
    m = go->model;

    gmxvcomm_loadxref(go, m->fnpdb);

    go->rhis_dx = m->rhis_dx;
    go->rhis_max = m->rhis_max;
    go->rhis_n = (int) ( go->rhis_max / go->rhis_dx + 0.5 );
    snew(go->rhis, go->rhis_n);

    go->wl = wl_openf(m->rmsdmin, m->rmsdmax, m->rmsddel,
        m->wl_lnf0, m->wl_flatness, m->wl_frac, m->invt_c, 0);
    fprintf(stderr, "parallel %d, domain-decomposition %d\n", PAR(cr), DOMAINDECOMP(cr));
  }

  if ( PAR(cr) ) {
    gmx_bcast(sizeof(gmxgomodel_t), go->model, cr);
  }

  go->kT = BOLTZ * tp;

  gmxgo_inithmc(go, mtop);

  return go;
}



static void gmxgo_close(gmxgo_t *go)
{
  if ( MASTER(go->cr) ) {
    sfree(go->rhis);
    wl_close(go->wl);
  }

  sfree(go->index);
  sfree(go->mass);
  if ( go->gvc != NULL ) {
    gmxvcomm_close(go->gvc);
  }
  sfree(go->x);
  sfree(go->f);
  sfree(go->xref);
  sfree(go->xf);
  sfree(go->x1);
  sfree(go->x2);
  sfree(go->stx);
  sfree(go->stv);
  sfree(go->stf);
  sfree(go);
}



static int gmxgo_saverhis(gmxgo_t *go, const char *fn)
{
  int i, imax;
  FILE *fp;

  for ( i = go->rhis_n - 1; i >= 0; i-- ) {
    if ( go->rhis[i] > 0 )
      break;
  }
  if ( i < 0 ) return -1;
  imax = i + 1;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  for ( i = 0; i < imax; i++ ) {
    fprintf(fp, "%g %g\n", go->rhis_dx * (i + 0.5), go->rhis[i]);
  }

  fprintf(stderr, "saving histogram file %s\n", fn);
  fclose(fp);
  return 0;
}



/* gather the coordinates `x` to `xg`,
 * which should `state->x` (local state)
 * only need for domain decomposition */
static int gmxgo_gatherx(gmxgo_t *go, rvec *x, double (*xg)[3], int doid)
{
  gmxvcomm_t *gvc = go->gvc;
  t_commrec *cr = go->cr;
  gmx_domdec_t *dd = cr->dd;
  int i, j, il, iw, lcnt, id, d;

  if ( !DOMAINDECOMP(cr) ) {
    /* direct map indices and return */
    for ( id = 0; id < go->n; id++ ) {
      for ( d = 0; d < 3; d++ ) {
        i = go->index[id];
        xg[id][d] = x[i][d];
      }
    }
    return 0;
  }

  if ( doid ) {
    /* we can save this call */
    gmxvcomm_gatherid(gvc);
  }

  /* fill the buffer `gvc->lx` with the local coordinates */
  for ( i = 0; i < gvc->lcnt; i++ ) {
    il = gvc->la[i];
    for ( d = 0; d < 3; d++ )
      gvc->lx[i][d] = x[il][d];
  }

  /* gather `lx` to `lx_m` in the master node */
  gmxvcomm_gatherv(gvc);

  /* extract the entire coordinates */
  if ( MASTER(cr) ) {
    /* clear the vector */
    for ( id = 0; id < go->n; id++ )
      vzero(xg[id]);

    for ( iw = 0, i = 0; i < dd->nnodes; i++ ) {
      if ( (lcnt = gvc->lcnt_m[i]) > 0 ) {
        for ( j = 0; j < lcnt; j++ ) {
          id = gvc->lwho_m[iw + j];
          for ( d = 0; d < 3; d++ )
            xg[id][d] = gvc->lx_m[iw + j][d];
          //printf("id %d: %g %g %g, from node %d\n",
          //    id, xg[id][0], xg[id][1], xg[id][2], i);
        }
        iw += lcnt;
      }
    }
  }

  return 0;
}



/* distribute the force additively */
static int gmxgo_scatterf(gmxgo_t *go, rvec *f, int doid)
{
  gmxvcomm_t *gvc = go->gvc;
  t_commrec *cr = go->cr;
  gmx_domdec_t *dd = cr->dd;
  int i, j, il, iw, lcnt, id, d;

  if ( !DOMAINDECOMP(cr) ) {
    /* direct map indices and return */
    for ( id = 0; id < go->n; id++ ) {
      for ( d = 0; d < 3; d++ ) {
        i = go->index[id];
        f[i][d] = (real) (f[i][d] + go->f[id][d]);
      }
    }
    return 0;
  }

  if ( doid ) {
    gmxvcomm_gatherid(gvc);
  }

  /* apply the force on the master node */
  if ( MASTER(cr) ) {
    for ( iw = 0, i = 0; i < dd->nnodes; i++ ) {
      if ( (lcnt = gvc->lcnt_m[i]) > 0 ) {
        for ( j = 0; j < lcnt; j++ ) {
          id = gvc->lwho_m[iw + j];
          vcopy(gvc->lx_m[iw + j], go->f[id]);
        }
        iw += lcnt;
      }
    }
  }

  /* gather `lx` to `lx_m` in the master node */
  gmxvcomm_scatterv(gvc);

  /* apply the force from the buffer `gvc->lx` */
  for ( i = 0; i < gvc->lcnt; i++ ) {
    il = gvc->la[i];
    /* apply the force additively */
    //fprintf(stderr, "i %d, %d: f (%g, %g, %g) += (%g, %g, %g)\n", gvc->lwho[i], gvc->ga[i], f[il][0], f[il][1], f[il][2], gvc->lx[i][0], gvc->lx[i][1], gvc->lx[i][2]);
    for ( d = 0; d < 3; d++ )
      f[il][d] = (real) ( f[il][d] + gvc->lx[i][d] );
  }

  return 0;
}



/* dump basic information */
__inline static void gmxgo_dump(gmxgo_t *go,
    double (*x)[3], int ePBC, matrix box, gmx_int64_t step)
{
  int i, c1, c;
  double dx[3], dr = 0, f = 0, dev = 0;
  t_pbc pbc[1];

  set_pbc(pbc, ePBC, box);

  fprintf(stderr, "step %s:\n", gmx_step_str(step, go->sbuf));
  for ( i = 0; i < go->n; i++ ) {
    c = gmxvcomm_where(go->gvc, i);
    if ( i > 0 ) {
      pbc_dx_d(pbc, x[i], x[i - 1], dx);
      c1 = gmxvcomm_where(go->gvc, i - 1);
      /* compute the consecutive peptide bond lengths */
      dr = vnorm(dx);
      f = vnorm(go->f[i]);
      dev = vdist(go->x[i], go->xf[i]);
    }
    fprintf(stderr, " %d[node %d]: f %g, dev %g, bond %g\n",
        i, c, f, dev, dr);
  }
  fprintf(stderr, "\n");
}



/* shift the coordinates to make the molecule whole
 * only the master needs to call this */
static int gmxgo_shift(gmxgo_t *go,
    double (*xin)[3], double (*xout)[3],
    int ePBC, matrix box, gmx_int64_t step)
{
  int i;
  double dx[3], dr;
  t_pbc pbc[1];

  set_pbc(pbc, ePBC, box);

  vcopy(xout[0], xin[0]);
  for ( i = 1; i < go->n; i++ ) {
    pbc_dx_d(pbc, xin[i], xin[i - 1], dx);
    dr = vnorm(dx); /* should be around 0.39 nm = 3.9 A */
    if ( dr > 0.45 || dr < 0.34 ) {
      fprintf(stderr, "step %s: %d-%d, dr %g, (%g, %g, %g); (%g, %g, %g) - (%g, %g, %g)\n",
          gmx_step_str(step, go->sbuf), i - 1, i, dr, dx[0], dx[1], dx[2],
          xin[i-1][0], xin[i-1][1], xin[i-1][2], xin[i][0], xin[i][1], xin[i][2]);
      gmxgo_dump(go, xin, ePBC, box, step);
      gmxgo_saverhis(go, go->model->fnrhis);
      wl_save(go->wl, go->model->fnvrmsd);
      return -1;
    }
    vadd(xout[i], xout[i - 1], dx);
  }
  return 0;
}



/* compute force */
static int gmxgo_calcf(gmxgo_t *go, double rmsd)
{
  int i;
  double dx[D], dvdx;
  gmxgomodel_t *m = go->model;

  dvdx = go->kT * wl_getdvf(go->wl, rmsd,
      m->mflmin, m->mflmax, m->mfhmin, m->mfhmax);
  if ( dvdx < m->mfmin ) {
    dvdx = m->mfmin;
  } else if ( dvdx > m->mfmax ) {
    dvdx = m->mfmax;
  }
  go->dvdx = dvdx;
  //fprintf(stderr, "rmsd %gA, dvdx %g\n", rmsd*10, go->dvdx);

  //dvdx  = go->kT * 2.0 * (rmsd - 0.5);
  dvdx /= rmsd * go->masstot;

  for ( i = 0; i < go->n; i++ )
    vzero( go->f[i] );

  for ( i = 0; i < go->n; i++ ) {
    vdiff(dx, go->x[i], go->xf[i]);
    vsinc(go->f[i], dx, -go->dvdx * go->mass[i]);
  }

  return 0;
}



/* ePBC == fr->ePBC */
static int gmxgo_rmsd_force(gmxgo_t *go, rvec *x, int doid, rvec *f,
    int ePBC, matrix box, gmx_int64_t step)
{
  double rmsd;
  int err = 0;
  gmxgomodel_t *m = go->model;
  t_commrec *cr = go->cr;

  /* collect `x` from different nodes to `go->x1` on the master */
  gmxgo_gatherx(go, x, go->x1, doid);

  if ( MASTER(cr) ) {
    /* make the coordinates whole */
    err = gmxgo_shift(go, go->x1, go->x, ePBC, box, step);
    if ( err != 0 ) goto BCAST_ERR;

    //for ( i = 0; i < go->n; i++ )
    //  printf("%d: %g %g %g | %g %g %g\n", i, go->x[i][0], go->x[i][1], go->x[i][2], go->xref[i][0], go->xref[i][1], go->xref[i][2]);

    rmsd = vrmsd(go->xref, go->xf, go->x, go->mass, go->n,
        0, NULL, NULL);
    wl_addf(go->wl, rmsd);
    wl_updatelnf(go->wl);

    /* update the histogram */
    if ( rmsd < go->rhis_max ) {
      go->rhis[ (int) ( rmsd / go->rhis_dx ) ] += 1;
    }

    if ( step > 0 && step % m->nstchat == 0 ) {
      fprintf(stderr, "\rstep %s: rmsd %gA, dvdx %g, flatness %g%%, lnf %g\n",
          gmx_step_str(step, go->sbuf), rmsd * 10, go->dvdx,
          wl_getflatness(go->wl) * 100, go->wl->lnf);
    }

    if ( step > 0 && step % m->nstrep == 0 ) {
      gmxgo_saverhis(go, m->fnrhis);
      wl_save(go->wl, m->fnvrmsd);
    }

    /* compute the force from the RMSD bias */
    gmxgo_calcf(go, rmsd);
  }

BCAST_ERR:
  if ( PAR(cr) ) {
    gmx_bcast(sizeof(err), &err, cr);
  }
  if ( err == 0 ) {
    gmxgo_scatterf(go, f, 0);
  }

  return err;
}



/* push the current x, v and f
 * we only save home atoms */
static int gmxgo_hmcpush(gmxgo_t *go,
    t_state *state, rvec *f)
{
  int i, n;
  t_commrec *cr = go->cr;

  /* determine the number of home atoms */
  n = PAR(cr) ? cr->dd->nat_home : state->natoms;

  if ( n > go->stncap ) {
    go->stncap = n;
    srenew(go->stx, n);
    srenew(go->stv, n);
    srenew(go->stf, n);
  }

  go->stn = n;
  memcpy(go->stx, state->x, sizeof(rvec) * n);
  memcpy(go->stv, state->v, sizeof(rvec) * n);
  memcpy(go->stf, f,        sizeof(rvec) * n);
  return 0;
}



/* pop x, v, f
 * assume the number of home atoms are not changed */
static void gmxgo_hmcpop(gmxgo_t *go,
    t_state *state, rvec *f, int reversev)
{
  int i, n = go->stn, nn, d;
  t_commrec *cr = go->cr;

  /* determine the number of home atoms */
  nn = PAR(cr) ? cr->dd->nat_home : state->natoms;
  if ( nn != n ) {
    fprintf(stderr, "nn %d vs. n %d\n", nn, n);
    exit(1);
  }

  if ( reversev ) {
    for ( i = 0; i < n; i++ ) {
      for ( d = 0; d < 3; d++ ) {
        go->stv[i][d] = -go->stv[i][d];
      }
    }
  }
  memcpy(state->x, go->stx, sizeof(rvec) * n);
  memcpy(state->v, go->stv, sizeof(rvec) * n);
  memcpy(f,        go->stf, sizeof(rvec) * n);
}



static double gmxgo_raw_rmsd(gmxgo_t *go,
    double (*xf)[3], double (*x)[3])
{
  int i;
  double dr2 = 0;

  for ( i = 0; i < go->n; i++ ) {
    dr2 += vdist2(xf[i], x[i]) * go->mass[i];
  }
  dr2 /= go->masstot;
  return sqrt( dr2 );
}



/* select */
static int gmxgo_hmcselect(gmxgo_t *go,
    t_state *state, rvec *f, int reversev,
    int ePBC, matrix box, gmx_int64_t step)
{
  double rmsd1, rmsd2, v1, v2;
  int err = 0, acc = 1;
  gmxgomodel_t *m = go->model;
  t_commrec *cr = go->cr;

  /* collect `state->x` from different nodes to `go->x1` on the master */
  gmxgo_gatherx(go, state->x, go->x1, 1);

  if ( MASTER(cr) ) {
    /* make the coordinates whole
     * TODO: possible shift due to PBC */
    err = gmxgo_shift(go, go->x1, go->x2, ePBC, box, step);
    if ( err != 0 ) goto BCAST_ERR;

    rmsd2 = vrmsd(go->xref, NULL, go->x2, go->mass, go->n,
        0, NULL, NULL);

    rmsd1 = gmxgo_raw_rmsd(go, go->xf, go->x2);

    v1 = wl_getvf(go->wl, rmsd1);
    v2 = wl_getvf(go->wl, rmsd2);

    acc = 1;
    if ( v2 > v1 ) {
      double r = rand01();
      if ( r > exp(-(v2 - v1)) ) {
        acc = 0;
        fprintf(stderr, "HMC rejecting a state, rmsd %g -> %g\n", rmsd1, rmsd2);
      }
    }
  }

BCAST_ERR:
  if ( PAR(cr) ) {
    gmx_bcast(sizeof(err), &err, cr);
    gmx_bcast(sizeof(acc), &acc, cr);
  }
  if ( !err && !acc ) {
    gmxgo_hmcpop(go, state, f, reversev);
  }

  return acc;
}



#endif /* GMXGO_H__ */
