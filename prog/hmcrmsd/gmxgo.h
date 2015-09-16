#ifndef GMXGO_H__
#define GMXGO_H__



/* Go model */
#define D 3
#include "mat.h"
#include "wl.h"
#include "hist.h"
#include "mtrand.h"
#include "gmxvcomm.h"
#include "gmxgomodel.h"



typedef struct {
  int n; /* number of atoms */
  int *index; /* index[0..n-1]: global indices of the Go atoms */
  double *mass; /* mass[0..n-1] mass */
  double masstot;
  double (*x)[D]; /* x[0..n-1] */
  double (*f)[D]; /* f[0..n-1] bias force */
  double (*xref)[D];
  double (*xf)[D]; /* fit structure */
  double (*x1)[D]; /* x1[0..n-1] */
  double (*x2)[D]; /* x2[0..n-1] */

  int nres; /* number of residues */
  int *resid;
  char **resnm;
  char **atnm;

  gmxvcomm_t *gvc;
  gmxgomodel_t model[1];
  double kT;
  wl_t *wl;
  hist_t *rhis; /* histogram */
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
  double hmcrej, hmctot;

  /* continue from the previous run */
  int isctn;

  FILE *fplog;
  /* string buffer to print step */
  char sbuf[STEPSTRSIZE];
} gmxgo_t;



/* block size used for dynamic allocation */
const int gmxgo_blksz = 32;



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

  /* set the number of residues */
  go->nres = mt->atoms.nres;

  /* find the molblock correpsonding to the moltype */
  id = 0;
  for ( ib = 0; ib < mtop->nmolblock; ib++ ) {
    mb = mtop->molblock + ib;
    if ( mb->type == itp ) break;
    id += mb->nmol * mb->natoms_mol;
  }

  /* search CA atoms from the moltype */
  go->n = 0;
  go->masstot = 0;
  ncap = gmxgo_blksz;
  snew(go->index, ncap);
  snew(go->mass,  ncap);
  snew(go->resid, ncap);
  snew(go->resnm, ncap);
  snew(go->atnm,  ncap);

  for ( ia = 0; ia < mt->atoms.nr; ia++ ) {
    char atnm[8];
    int resind = mt->atoms.atom[ia].resind;
    int sel = 0, seltype = go->model->seltype;

    strcpy(atnm, mt->atoms.atomname[ia][0]);
    strstrip(atnm);

    if ( seltype == SEL_CA ) {
      if ( strcmp(atnm, "CA") == 0 ) {
        sel = 1;
      }
    } else if ( seltype == SEL_HEAVY ) {
      if ( atnm[0] != 'H' && atnm[0] != 'D' ) {
        sel = 1;
      }
    } else if ( seltype == SEL_ALL ) {
      sel = 1;
    } else if ( seltype == SEL_CAENDTOEND ) {
      if ( strcmp(atnm, "CA") == 0
        &&  (resind == 0 || resind == go->nres - 1) ) {
        sel = 1;
      }
    } else {
      fprintf(stderr, "unknown seltype %d\n", seltype);
      return -1;
    }

    if ( sel ) {
      /* reallocate the memory */
      if ( go->n >= ncap ) {
        ncap += gmxgo_blksz;
        srenew(go->index, ncap);
        srenew(go->mass,  ncap);
        srenew(go->resid, ncap);
        srenew(go->resnm, ncap);
        srenew(go->atnm,  ncap);
      }

      go->index[go->n] = id + ia;
      go->mass[go->n] = mt->atoms.atom[ia].m;
      go->masstot += go->mass[go->n];
      /* `resind` is the zero-based index
       * `resinfo[].nr` is the residue index
       * that is supposed to match that in PDB */
      go->resid[go->n] = mt->atoms.resinfo[resind].nr;
      go->resnm[go->n] = mt->atoms.resinfo[resind].name[0];
      go->atnm[go->n] = mt->atoms.atomname[ia][0];
      go->n += 1;
    }
  }

  fprintf(stderr, "%d spatoms, total mass %g:\n",
      go->n, go->masstot);
  for ( i = 0; i < go->n; i++ ) {
    fprintf(stderr, " %4d %-4s %4d %-4s %8.3f\n",
        go->index[i], go->atnm[i], go->resid[i], go->resnm[i],
        go->mass[i]);
  }

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
static int gmxgo_loadxref(gmxgo_t *go, const char *fnpdb)
{
  FILE *fp;
  char s[256];

  int i, id, err;
  int natm, ncap, *resid, *used, *matched;
  float (*xref)[D], x[D];
  char (*atnm)[8], (*resnm)[8], sresid[8];

  natm = 0;
  ncap = gmxgo_blksz;
  snew(resid, ncap);
  snew(xref,  ncap);
  snew(atnm,  ncap);
  snew(resnm, ncap);

  if ( (fp = fopen(fnpdb, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fnpdb);
    return -1;
  }

  /* 1. load raw PDB information */
  while ( fgets(s, sizeof s, fp) ) {
    /* a model is finished, we implicitly include ENDMDL */
    if ( strncmp(s, "TER", 3) == 0 || strncmp(s, "END", 3) == 0 )
      break;

    /* PDB format specification
     * http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
     * */
    if ( strncmp(s, "ATOM  ", 6) != 0 )
      continue;

    /* alternate location */
    if ( s[16] != ' ' && s[16] != 'A' )
      continue;

    /* reallocate memory if possible */
    if ( natm >= ncap ) {
      ncap += gmxgo_blksz;
      srenew(resid, ncap);
      srenew(xref,  ncap);
      srenew(atnm,  ncap);
      srenew(resnm, ncap);
    }

    /* copy the atom name */
    strncpy(atnm[natm], s + 12, 4);
    atnm[natm][4] = '\0';
    strstrip(atnm[natm]);

    /* copy the residue name */
    strncpy(resnm[natm], s + 17, 3);
    resnm[natm][3] = '\0';
    strstrip(resnm[natm]);

    /* copy the residue index */
    strncpy(sresid, s + 22, 4);
    sresid[4] = '\0';
    resid[natm] = atoi( strstrip(sresid) );

    /* scan the coordinates */
    sscanf(s + 30 , "%f%f%f", x, x + 1, x + 2);

    /* convert the unit from angstroms to nanometers */
    xref[natm][0] = x[0] / 10;
    xref[natm][1] = x[1] / 10;
    xref[natm][2] = x[2] / 10;

    /*
    fprintf(stderr, "%4d %-4s %4d %-4s %8.3f %8.3f %8.3f\n",
        natm, atnm[natm], resid[natm], resnm[natm],
        xref[natm][0], xref[natm][1], xref[natm][2]);
    */
    natm += 1;
  }

  /* 2. match raw PDB information to the GROMACS topology */
  err = 0;
  /* the array `used` tracks if each PDB entry is used */
  snew(used,  natm);
  for ( i = 0; i < natm; i++ ) used[i] = 0;

  snew(matched, go->n);
  for ( id = 0; id < go->n; id++ ) matched[id] = 0;

  for ( id = 0; id < go->n; id++ ) {
    /* try to find the matching atom from the PDB file */
    for ( i = 0; i < natm; i++ ) {
      if ( resid[i] == go->resid[id]
        && strcmp(resnm[i], go->resnm[id]) == 0
        && strcmp(atnm[i], go->atnm[id]) == 0 ) {
        break;
      }
    }

    /* if matching fails */
    if ( i >= natm ) {
      if ( err <= 1
        && ( go->atnm[id][0] == 'H' || go->atnm[id][0] == 'D' ) ) {
        err = 1;
      } else {
        fprintf(stderr, "cannot find matching atom for "
            "id %d, atnm %s, resnm %s, resid %d\n",
            id, go->atnm[id], go->resnm[id], go->resid[id]);
        err = 2;
        break;
      }
    } else {
      matched[id] = 1;
    }

    /* copy the reference coordinates */
    go->xref[id][0] = xref[i][0];
    go->xref[id][1] = xref[i][1];
    go->xref[id][2] = xref[i][2];

    /* make this PDB entry is used, so we avoid
     * this in subsequent fuzzy matching */
    used[i] = 1;
  }

  /* try fuzzy matching for hydrogen or deuterium atoms */
  if ( err == 1 ) {
    int len;

    err = 0; /* clear the error flag */
    for ( id = 0; id < go->n; id++ ) {
      if ( matched[id] ) continue;

      len = strlen(go->atnm[id]) - 1;

      /* we do fuzzy matching only if the last character
       * of the atom name is a number */
      if ( !isdigit(go->atnm[id][len]) ) {
        err = 2;
      } else {
        for ( i = 0; i < natm; i++ ) {
          if ( resid[i] == go->resid[id]
            && strcmp(resnm[i], go->resnm[id]) == 0
            && strncmp(atnm[i], go->atnm[id], len) == 0
            && !used[i] ) {
            fprintf(stderr, "fuzzy matching atom %4d, "
                "%4s (PDB) -> %4s (GROMACS, residue %s%d)\n",
                id, atnm[i], go->atnm[id],
                go->resnm[i], go->resid[id]);
            break;
          }
        }

        if ( i >= natm ) {
          err = 2;
        } else {
          /* copy the reference coordinates */
          go->xref[id][0] = xref[i][0];
          go->xref[id][1] = xref[i][1];
          go->xref[id][2] = xref[i][2];

          /* mark this PDB entry as used, so we can avoid
           * this in subsequent fuzzy matching */
          used[i] = 1;
          matched[id] = 1;
        }
      }

      if ( err ) {
        fprintf(stderr, "cannot fuzzy match atom for "
            "id %4d, atnm %4s, resnm %4s, resid %4d\n",
            id, go->atnm[id], go->resnm[id], go->resid[id]);
        break;
      }
    }
  }

  free(resid);
  free(used);
  free(matched);
  free(xref);
  free(atnm);
  free(resnm);

  return err;
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
    n = (int) ((n + 3.0 * sqrt(n)) / cr->nnodes);
  }

  go->stncap = n;
  go->stn = 0;
  snew(go->stx, n);
  snew(go->stv, n);
  snew(go->stf, n);
  go->hmcrej = 0;
  go->hmctot = 0;
}



/* tp: the temperature
 * ctn: continue from the previous simulation */
static gmxgo_t *gmxgo_open(gmx_mtop_t *mtop, t_commrec *cr,
    double tp, const char *fncfg, int ctn)
{
  gmxgo_t *go;

  snew(go, 1);

  go->cr = cr;

  /* load input parameters from the configuration file */
  if ( MASTER(cr) ) {
    gmxgomodel_t *m = go->model;

    gmxgomodel_default(m);
    gmxgomodel_load(m, fncfg);
  }

  /* broadcast the parameters */
  if ( PAR(cr) ) {
    gmx_bcast(sizeof(gmxgomodel_t), go->model, cr);
  }

  go->kT = BOLTZ * tp;


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

    if ( gmxgo_loadxref(go, m->fnpdb) != 0 ) {
      exit(1);
    }

    go->wl = wl_openf(m->rmsdmin, m->rmsdmax, m->rmsddel,
        m->wl_lnf0, m->wl_flatness, m->wl_frac, m->invt_c, 0);

    go->rhis = hist_open(0, m->rhis_max, m->rhis_dx);

    /* continue from the last run */
    go->isctn = ctn;
    if ( ctn ) {
      wl_load(go->wl, m->fnvrmsd);
      hist_load(go->rhis, m->fnrhis);
    }
    fprintf(stderr, "parallel %d, domain-decomposition %d\n",
        PAR(cr), DOMAINDECOMP(cr));
  }

  gmxgo_inithmc(go, mtop);

  return go;
}



static void gmxgo_close(gmxgo_t *go)
{
  if ( MASTER(go->cr) ) {
    wl_close(go->wl);
    hist_close(go->rhis);
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



/* gather the coordinates `x` to `xg`,
 * which should `state->x` (local state)
 * only need for domain decomposition */
static int gmxgo_gatherx(gmxgo_t *go, rvec *x,
    double (*xg)[D], int doid)
{
  gmxvcomm_t *gvc = go->gvc;
  t_commrec *cr = go->cr;
  gmx_domdec_t *dd = cr->dd;
  int i, j, il, iw, lcnt, id, d;

  if ( !DOMAINDECOMP(cr) ) {
    /* direct map indices and return */
    for ( id = 0; id < go->n; id++ ) {
      i = go->index[id]; /* global index */
      for ( d = 0; d < D; d++ ) {
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
    for ( d = 0; d < D; d++ )
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
        for ( j = 0; j < lcnt; j++, iw++ ) {
          id = gvc->lwho_m[iw]; /* special-atom index */
          for ( d = 0; d < D; d++ )
            xg[id][d] = gvc->lx_m[iw][d];
          //printf("id %d: %g %g %g, from node %d\n",
          //    id, xg[id][0], xg[id][1], xg[id][2], i);
        }
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
      i = go->index[id]; /* global index */
      for ( d = 0; d < D; d++ ) {
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
        for ( j = 0; j < lcnt; j++, iw++ ) {
          id = gvc->lwho_m[iw]; /* special-atom index */
          vcopy(gvc->lx_m[iw], go->f[id]);
        }
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
    for ( d = 0; d < D; d++ )
      f[il][d] = (real) ( f[il][d] + gvc->lx[i][d] );
  }

  return 0;
}



/* distance between two alpha carbon atoms */
#define CACABOND 0.382
/* maximal deviation */
#define CACABOND_DEV 0.06



/* dump basic information */
__inline static void gmxgo_dump(gmxgo_t *go,
    double (*x)[D], int ePBC, matrix box, gmx_int64_t step)
{
  int i, c1, c;
  double dx[D], dr = 0, f = 0, dev = 0, rmsd;
  t_pbc pbc[1];

  set_pbc(pbc, ePBC, box);

  rmsd = vrmsd(go->xref, go->xf, go->x, go->mass, go->n,
      0, NULL, NULL);
  fprintf(stderr, "step %s, RMSD %gA:\n",
      gmx_step_str(step, go->sbuf), rmsd * 10);

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
    fprintf(stderr, " %d[node %d]: f %g, dev %gA, bond %gA (%+gA)\n",
        i, c, f, dev * 10, dr * 10, (dr - CACABOND) * 10);
  }
  fprintf(stderr, "\n");
}



/* shift the coordinates to make the molecule whole
 * only the master needs to call this */
static int gmxgo_shift(gmxgo_t *go,
    double (*xin)[D], double (*xout)[D],
    int ePBC, matrix box, gmx_int64_t step)
{
  int i, id = -1, err = 0;
  double dx[D], dr, dev, devmax = 0;
  t_pbc pbc[1];

  set_pbc(pbc, ePBC, box);

  vcopy(xout[0], xin[0]);
  for ( i = 1; i < go->n; i++ ) {
    pbc_dx_d(pbc, xin[i], xin[i - 1], dx);
    if ( go->model->seltype == SEL_CA ) {
      dr = vnorm(dx);
      dev = fabs(dr - CACABOND);
      if ( dev > devmax ) {
        devmax = dev;
        id = i;
      }
    }
    vadd(xout[i], xout[i - 1], dx);
  }

  err = (devmax > CACABOND_DEV);
  if ( err ) {
    pbc_dx_d(pbc, xin[id], xin[id - 1], dx);
    dr = vnorm(dx);
    fprintf(stderr, "step %s: broken peptide bond %d-%d, dr %g, "
        "(%g, %g, %g); (%g, %g, %g) - (%g, %g, %g), dev %g\n",
        gmx_step_str(step, go->sbuf), id - 1, id, dr,
        dx[0], dx[1], dx[2],
        xin[id - 1][0], xin[id - 1][1], xin[id - 1][2],
        xin[id][0], xin[id][1], xin[id][2], devmax);
    gmxgo_dump(go, xin, ePBC, box, step);
    hist_save(go->rhis, go->model->fnrhis);
    wl_save(go->wl, go->model->fnvrmsd);
  }

  //fprintf(stderr, "step %d, max dev %gA\n", (int) step, devmax*10);
  return err;
}



/* project the force `f` along rmsd */
static double gmxgo_projectf(gmxgo_t *go, double (*f)[D],
    double rmsd, double kT, double (*x)[D], double (*xf)[D])
{
  int i, n = go->n;
  double dx[D], dot = 0;

  for ( i = 0; i < n; i++ ) {
    /* compute the gradient */
    vdiff(dx, x[i], xf[i]);
    dot += vdot(dx, f[i]);
  }

  /* because of the fixed length of the peptide bonds
   * the factor D * n may not be correct */
  return (dot / kT + D * n - 1) / rmsd;
}



/* compute the mean force `go->mf` of the current `rmsd`
 * does not apply to the actual force */
static int gmxgo_calcmf(gmxgo_t *go, double rmsd)
{
  int i;
  double dx[D], dvdx;
  gmxgomodel_t *m = go->model;

  if ( m->bias_mf ) {
    /* compute the gradient from the bias potential */
    dvdx = wl_getdvdx_mf(go->wl, rmsd, m->minh,
        m->mflmin, m->mflmax, m->mfhmin, m->mfhmax);
    //printf("%g, %g\n", rmsd, dvdx); getchar();
  } else {
    /* compute the gradient from the bias mean force */
    dvdx = wl_getdvdx_v(go->wl, rmsd,
        m->mflmin, m->mflmax, m->mfhmin, m->mfhmax);
  }
  dvdx *= go->kT;

  /* place the limits */
  if ( dvdx < m->mfmin ) {
    dvdx = m->mfmin;
  } else if ( dvdx > m->mfmax ) {
    dvdx = m->mfmax;
  }

  go->dvdx = dvdx;
  //fprintf(stderr, "rmsd %gA, dvdx %g\n", rmsd*10, go->dvdx);

  //dvdx  = go->kT * 2.0 * (rmsd - 0.5);
  dvdx /= rmsd * go->masstot;

  for ( i = 0; i < go->n; i++ ) {
    vzero( go->f[i] );
  }

  for ( i = 0; i < go->n; i++ ) {
    vdiff(dx, go->x[i], go->xf[i]);
    vsinc(go->f[i], dx, -go->dvdx * go->mass[i]);
  }

  return 0;
}



/* print a short message on screen */
static void gmxgo_chat(gmxgo_t *go, gmx_int64_t step, double rmsd)
{
  fprintf(stderr, "\rstep %s: rmsd %gA, dvdx %g, flatness %g%%, lnf %g, hmcrej %g\n",
      gmx_step_str(step, go->sbuf), rmsd * 10, go->dvdx,
      wl_getflatness(go->wl) * 100, go->wl->lnf, go->hmcrej);
}



/* write the log file */
static void gmxgo_log(gmxgo_t *go, gmx_int64_t step, double rmsd)
{
  if ( go->fplog == NULL ) {
    go->fplog = fopen(go->model->fnlog, (go->isctn ? "a" : "w"));
  }

  if ( go->fplog != NULL ) {
    fprintf(go->fplog, "%s %g\n",
        gmx_step_str(step, go->sbuf), rmsd);
  }
}



/* ePBC == fr->ePBC */
static int gmxgo_rmsd_force(gmxgo_t *go, t_state *state, int doid, rvec *f,
    int ePBC, matrix box, gmx_int64_t step)
{
  double rmsd = 0, mf;
  rvec *x = state->x;
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
    hist_add(go->rhis, rmsd);

    if ( (step > 0 && step % m->nstchat == 0) || m->debug ) {
      gmxgo_chat(go, step, rmsd);

      /* debugging code to test HMC rejection */
      if ( m->debug == 3 ) {
        int id;

        for ( id = 0; id < go->stn; id+= 1000 ) {
          fprintf(stderr, "step %d, id %d, x: %g, %g, %g\n", (int) step, id, x[id][0], x[id][1], x[id][2]);
          fprintf(stderr, "step %d, id %d, v: %g, %g, %g\n", (int) step, id, state->v[id][0], state->v[id][1], state->v[id][2]);
          fprintf(stderr, "step %d, id %d, f: %g, %g, %g\n", (int) step, id, f[id][0], f[id][1], f[id][2]);
        }
        //getchar();
      }

      if ( step > 0 && step % m->nstlog == 0 ) {
        gmxgo_log(go, step, rmsd);
      }
    }

    if ( step > 0 && step % m->nstrep == 0 ) {
      wl_save(go->wl, m->fnvrmsd);
      hist_save(go->rhis, m->fnrhis);
    }
  }

BCAST_ERR:
  if ( PAR(cr) ) {
    gmx_bcast(sizeof(err), &err, cr);
  }
  if ( err ) return err;

  /* compute the projection of force on RMSD
   * we compute the mean force even it is not needed */
  /* if ( m->bias_mf ) */
  {
    /* collect `f` from different nodes to `go->x1` on the master
     * `doid = 0` since we must have indices now */
    gmxgo_gatherx(go, f, go->x1, 0);

    if ( MASTER(cr) ) {
      /* `go->x1` is the force, project it along RMSD
       * note that the force must be the force without the bias */
      mf = gmxgo_projectf(go, go->x1, rmsd, go->kT, go->x, go->xf);
      //fprintf(stderr, "step %d: computing mean force, rmsd %g, mf %g\n", (int) step, rmsd, mf); getchar();
      wl_addforcef(go->wl, rmsd, mf);
    }
  }

  if ( MASTER(cr) ) {
    /* compute the force from the RMSD bias */
    gmxgo_calcmf(go, rmsd);
  }

  /* apply the mean force */
  gmxgo_scatterf(go, f, 0);

  return err;
}



/* copy n vectors */
static void gmxgo_vcopy(rvec *dest, rvec *src, int n)
{
  int i, d = 0;

  //memcpy(dest, src, sizeof(rvec) * n);
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < DIM; d++ ) {
      dest[i][d] = src[i][d];
    }
  }
}


/* push the current x and f
 * we only save home atoms
 * call this after do_force() to make sure
 * that the force `f` matches the position `x` */
static int gmxgo_hmcpushxf(gmxgo_t *go,
    t_state *state, rvec *f, gmx_int64_t step)
{
  int i, n;
  t_commrec *cr = go->cr;

  /* determine the number of home atoms */
  n = PAR(cr) ? cr->dd->nat_home : state->natoms;

  if ( n > go->stncap ) {
    if ( MASTER(cr) ) {
      fprintf(stderr, "hmcpushxf: step %s, expanding state variable %d -> %d\n",
          gmx_step_str(step, go->sbuf), go->stncap, n);
    }
    go->stncap = n;
    srenew(go->stx, n);
    srenew(go->stv, n);
    srenew(go->stf, n);
  }

#if 0
  if ( go->model->debug == 3 && step % 300 - 1 == 3 ) {
    int id;

    /* see if the force matches the just popped value */
    for ( id = 0; id < n; id+= 1000 ) {
      fprintf(stderr, "step %d, id %d, f: %g(%g), %g(%g), %g(%g)\n",
          (int) step, id,
          f[id][0], go->stf[id][0],
          f[id][1], go->stf[id][1],
          f[id][2], go->stf[id][2]);
    }
    getchar();
  }
#endif

  go->stn = n;
  gmxgo_vcopy(go->stx, state->x, n);
  gmxgo_vcopy(go->stf, f,        n);

  return 0;
}



/* push the current v
 * we only save home atoms */
static int gmxgo_hmcpushv(gmxgo_t *go,
    t_state *state, gmx_int64_t step)
{
  int i, n;
  t_commrec *cr = go->cr;

  /* determine the number of home atoms */
  n = PAR(cr) ? cr->dd->nat_home : state->natoms;

  if ( n > go->stncap ) {
    if ( MASTER(cr) ) {
      fprintf(stderr, "hmcpushv: step %s, expanding state variable %d -> %d\n",
          gmx_step_str(step, go->sbuf), go->stncap, n);
    }
    go->stncap = n;
    srenew(go->stx, n);
    srenew(go->stv, n);
    srenew(go->stf, n);
  }

  go->stn = n;
  gmxgo_vcopy(go->stv, state->v, n);
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
      for ( d = 0; d < D; d++ ) {
        go->stv[i][d] = -go->stv[i][d];
      }
    }
  }

  gmxgo_vcopy(state->x, go->stx, n);
  gmxgo_vcopy(state->v, go->stv, n);
  gmxgo_vcopy(f,        go->stf, n);
}



static double gmxgo_raw_rmsd(gmxgo_t *go,
    double (*xf)[D], double (*x)[D])
{
  int i;
  double dr2, s = 0;

  for ( i = 0; i < go->n; i++ ) {
    dr2 = vdist2(xf[i], x[i]);
    s += dr2 * go->mass[i];
  }
  return sqrt( s / go->masstot );
}



/* select */
static int gmxgo_hmcselect(gmxgo_t *go,
    t_state *state, rvec *f, int reversev,
    int ePBC, matrix box, gmx_int64_t step)
{
  double rmsd1, rmsd2, delv, dr0;
  int err = 0, acc = 1;
  gmxgomodel_t *m = go->model;
  t_commrec *cr = go->cr;

  /* collect `state->x` from different nodes to `go->x1` on the master */
  gmxgo_gatherx(go, state->x, go->x1, 1);

  if ( MASTER(cr) ) {
    /* remove the periodic boundary condition
     * to make the coordinates whole */
    err = gmxgo_shift(go, go->x1, go->x2, ePBC, box, step);
    if ( err != 0 ) goto BCAST_ERR;

    /* transform `xref` to `x1` to fit `x2`
     * compute the RMSD between `x1` and `x2` */
    rmsd2 = vrmsd(go->xref, go->x1, go->x2, go->mass, go->n,
        0, NULL, NULL);

    /* comupte the RMSD between `xf` and `x2` */
    rmsd1 = gmxgo_raw_rmsd(go, go->xf, go->x2);

    /* the two fitted structures are supposed to be close
     * if not, one of them might have been wrapped unintendedly */
    dr0 = vdist(go->xf[0], go->x1[0]);
    if ( dr0 > 0.1 ) {
      /* the threshold 0.1nm or 1A is generous, can be smaller */
      fprintf(stderr, "The first special atoms of the two "
          "fitted structures differ by %gA\n", dr0 * 10);
      goto BCAST_ERR;
    }

    if ( m->bias_mf ) {
      delv = wl_getdelv_mf(go->wl, rmsd1, rmsd2, m->minh,
          m->mflmin, m->mflmax, m->mfhmin, m->mfhmax);
    } else {
      delv = wl_getdelv_v(go->wl, rmsd1, rmsd2);
    }

    acc = 1;
    if ( delv > 0 ) {
      double r = rand01();
      acc = ( r < exp(-delv) );
    }

    /* special flag to test HMC rejection */
    if ( m->debug == 3 ) {
      acc = ( step % 300 != 3 );
    }

    if ( !acc ) {
      go->hmcrej += 1;
      fprintf(stderr, "step %s: HMC rejecting a state, rmsd %g -> %g, delv %g\n",
          gmx_step_str(step, go->sbuf), rmsd1, rmsd2, delv);
    }
    go->hmctot += 1;
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
