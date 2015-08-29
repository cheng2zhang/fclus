#ifndef GMXGO_H__
#define GMXGO_H__



/* Go model */
#define D 3
#include "mat.h"
#include "gmxvcomm.h"




typedef struct {
  int n; /* number of atoms */
  int *index; /* index[0..n-1]: global indices of the Go atoms */
  double (*x)[3]; /* x[0..n-1] */
  double (*f)[3]; /* f[0..n-1] bias force */
  double (*xref)[3];
  double (*xf)[3]; /* fit structure */
  double (*x1)[3]; /* x1[0..n-1] */
  gmxvcomm_t *gvc;
  double kT;
  double *rhis; /* histogram */
  double rhis_dx, rhis_max; /* spacing and maximal */
  int rhis_n;
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
  ncap = BLKSZ;
  snew(go->index, ncap);
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
      }

      go->index[go->n] = id + ia;
      go->n += 1;
    }
  }

  fprintf(stderr, "%d spatoms:", go->n);
  for ( i = 0; i < go->n; i++ ) {
    fprintf(stderr, " %d", go->index[i]);
  }
  fprintf(stderr, "\n");

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
    }
    gmx_bcast(sizeof(go->index[0]) * go->n, go->index, cr);
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



/* tp: the temperature */
static gmxgo_t *gmxgo_open(gmx_mtop_t *mtop, t_commrec *cr,
    double tp)
{
  gmxgo_t *go;

  snew(go, 1);

  gmxgo_build(go, mtop, cr);
  go->gvc = gmxvcomm_open(go->n, go->index, cr);

  snew(go->x, go->n);
  snew(go->f, go->n);
  snew(go->xref, go->n);
  snew(go->xf, go->n);
  snew(go->x1, go->n);

  go->kT = BOLTZ * tp;
  /* load the reference */
  if ( MASTER(cr) ) {
    gmxvcomm_loadxref(go, "1VII.pdb");

    go->rhis_dx = 0.002;
    go->rhis_max = 2.0; /* 20 angstrom */
    go->rhis_n = (int) ( go->rhis_max / go->rhis_dx + 0.5 );
    snew(go->rhis, go->rhis_n);
  }

  return go;
}



static void gmxgo_close(gmxgo_t *go)
{
  if ( MASTER(go->gvc->cr) ) {
    sfree(go->rhis);
  }

  sfree(go->index);
  gmxvcomm_close(go->gvc);
  sfree(go->x);
  sfree(go->f);
  sfree(go->xref);
  sfree(go->xf);
  sfree(go->x1);
  sfree(go);
}



/* gather the coordinates `x` to `xg`,
 * which should `state->x` (local state) */
static int gmxgo_gatherx(gmxgo_t *go, rvec *x, double (*xg)[3], int doid)
{
  gmxvcomm_t *gvc = go->gvc;
  t_commrec *cr = gvc->cr;
  gmx_domdec_t *dd = cr->dd;
  int i, j, il, iw, lcnt, id, d;

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



/* distribute the force */
static int gmxgo_scatterf(gmxgo_t *go, rvec *f, int doid)
{
  gmxvcomm_t *gvc = go->gvc;
  t_commrec *cr = gvc->cr;
  gmx_domdec_t *dd = cr->dd;
  int i, j, il, iw, lcnt, id, d;

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

  //if ( MASTER(cr) ) getchar();
  //gmx_barrier(gvc->cr);
  return 0;
}



/* shift the coordinates to make the molecule whole */
static void gmxgo_shift(gmxgo_t *go,
    double (*xin)[3], double (*xout)[3],
    int ePBC, matrix box)
{
  int i;
  double dx[3];
  t_pbc pbc[1];

  set_pbc(pbc, ePBC, box);

  vcopy(xout[0], xin[0]);
  for ( i = 1; i < go->n; i++ ) {
    pbc_dx_d(pbc, xin[i], xin[i - 1], dx);
    vadd(xout[i], xout[i - 1], dx);
  }
}



/* compute force */
static int gmxgo_calcf(gmxgo_t *go, double rmsd)
{
  int i;
  double dx[D], dvdx;

  //dvdx = go->kT * wl_getdvf(wl, rmsd, mfl, mfh);
  dvdx  = go->kT * 10.0 * (rmsd - 0.5);
  dvdx /= rmsd * go->n;

  for ( i = 0; i < go->n; i++ ) {
    vdiff(dx, go->x[i], go->xf[i]);
    vsinc(go->f[i], dx, -dvdx);
  }

  return 0;
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



/* ePBC == fr->ePBC */
static int gmxgo_rmsd(gmxgo_t *go, rvec *x, int doid, rvec *f,
    int ePBC, matrix box, gmx_int64_t step)
{
  double rmsd;
  int i;
  char sbuf[STEPSTRSIZE];

  /* collect `x` from different nodes to `go->x1` on the master */
  gmxgo_gatherx(go, x, go->x1, doid);

  /* make the coordinates whole */
  gmxgo_shift(go, go->x1, go->x, ePBC, box);

  if ( MASTER(go->gvc->cr) ) {
    //for ( i = 0; i < go->n; i++ )
    //  printf("%d: %g %g %g | %g %g %g\n", i, go->x[i][0], go->x[i][1], go->x[i][2], go->xref[i][0], go->xref[i][1], go->xref[i][2]);
    rmsd = vrmsd(go->xref, go->xf, go->x, NULL, go->n,
        0, NULL, NULL);
    if ( rmsd < go->rhis_max ) {
      go->rhis[ (int) ( rmsd / go->rhis_dx ) ] += 1;
    }
    if ( step > 0 && step % 10000 == 0 ) {
      gmxgo_saverhis(go, "rhis.dat");
    }

    fprintf(stderr, "step %s: rmsd %g A\n", gmx_step_str(step, sbuf), rmsd * 10);
    gmxgo_calcf(go, rmsd);
  }

  gmxgo_scatterf(go, f, 0);

  return 0;
}



#endif /* GMXGO_H__ */
