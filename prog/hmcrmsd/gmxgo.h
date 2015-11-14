#ifndef GMXGO_H__
#define GMXGO_H__



/* Go model */
#include "util.h"
#define D 3
#include "mat.h"
#include "wl.h"
#include "hist.h"
#include <time.h>
#include "mtrand.h"
#include "gmxvcomm.h"
#include "gmxgocfg.h"



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
  double (*xwhole)[D];
  double (*xwholep)[D];
  double (*xrt)[D];
  double (*xrtp)[D];

  int nres; /* number of residues */
  int *resid; /* resid[0..n-1] is the PDB residue ID */
  char **resnm;
  char **atnm;

  gmxvcomm_t *gvc;
  gmxgocfg_t cfg[1];

  double kT;
  double *wr;
  wl_t *wl;

  hist_t *rhis; /* histogram */
  t_commrec *cr;
  double dvdx;

  double rmsd;

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
const int gmxgo_blksz = 256;



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

  /* find the molblock corresponding to the moltype */
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
  xnew(go->index, ncap);
  xnew(go->mass,  ncap);
  xnew(go->resid, ncap);
  xnew(go->resnm, ncap);
  xnew(go->atnm,  ncap);

  for ( ia = 0; ia < mt->atoms.nr; ia++ ) {
    char atnm[8];
    int resind = mt->atoms.atom[ia].resind;
    int sel = 0, rmsdgrp = go->cfg->rmsdgrp;

    strcpy(atnm, mt->atoms.atomname[ia][0]);
    strstrip(atnm);

    if ( rmsdgrp == RMSDGRP_CA ) {
      if ( strcmp(atnm, "CA") == 0 ) {
        sel = 1;
      }
    } else if ( rmsdgrp == RMSDGRP_HEAVY ) {
      if ( atnm[0] != 'H' && atnm[0] != 'D' ) {
        sel = 1;
      }
    } else if ( rmsdgrp == RMSDGRP_ALL ) {
      sel = 1;
    } else if ( rmsdgrp == RMSDGRP_CAENDTOEND ) {
      if ( strcmp(atnm, "CA") == 0
        &&  (resind == 0 || resind == go->nres - 1) ) {
        sel = 1;
      }
    } else {
      fprintf(stderr, "unknown rmsdgrp %d\n", rmsdgrp);
      return -1;
    }

    if ( sel ) {
      /* reallocate the memory */
      if ( go->n >= ncap ) {
        ncap += gmxgo_blksz;
        xrenew(go->index, ncap);
        xrenew(go->mass,  ncap);
        xrenew(go->resid, ncap);
        xrenew(go->resnm, ncap);
        xrenew(go->atnm,  ncap);
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



/* extract the coordinates of the special atoms
 * from the GROMACS topology */
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
      xnew(go->index, go->n);
      xnew(go->mass, go->n);
    }
    gmx_bcast(sizeof(go->index[0]) * go->n, go->index, cr);
    gmx_bcast(sizeof(go->mass[0]) * go->n, go->mass, cr);
    gmx_bcast(sizeof(go->masstot), &go->masstot, cr);
  }

  return 0;
}



/* load the reference coordinates from a PDB file */
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
  /* the size of the following arrays are increased in multiples
   * of gmxgo_blksz */
  xnew(resid, ncap);
  xnew(xref,  ncap);
  xnew(atnm,  ncap);
  xnew(resnm, ncap);

  if ( (fp = fopen(fnpdb, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fnpdb);
    return -1;
  }

  /* 1. load all atom information from the PDB file */
  while ( fgets(s, sizeof s, fp) ) {
    /* stop when the first model is finished,
     * we implicitly include ENDMDL */
    if ( strncmp(s, "TER", 3) == 0
      || strncmp(s, "END", 3) == 0 ) {
      break;
    }

    /* PDB format specification
     * http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
     * */
    if ( strncmp(s, "ATOM  ", 6) != 0 )
      continue;

    /* alternate location */
    if ( s[16] != ' ' && s[16] != 'A' )
      continue;

    /* reallocate memory if needed */
    if ( natm >= ncap ) {
      ncap += gmxgo_blksz;
      xrenew(resid, ncap);
      xrenew(xref,  ncap);
      xrenew(atnm,  ncap);
      xrenew(resnm, ncap);
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

  /* `used[i]` tracks if the `i`th PDB ATOM entry is used */
  xnew(used,  natm);
  for ( i = 0; i < natm; i++ ) used[i] = 0;

  /* `matched[id]` tracks if the `id`th needed atom has found
   * a correspondence in the PDB ATOM entries */
  xnew(matched, go->n);
  for ( id = 0; id < go->n; id++ ) matched[id] = 0;

  for ( id = 0; id < go->n; id++ ) {
    /* try to find the matching atom from the PDB file
     * for the `id`th needed atom */
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

    /* mark this PDB entry as used, so we can avoid
     * it in subsequent fuzzy matching */
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
       * of the atom name is a number, like HA2 */
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



/* build the target distribution function */
static double *gmxgo_mkwr(double xmin, double xmax, double dx,
    double exponent)
{
  int i, n;
  double x, *wr;

  n = (int) ( (xmax - xmin) / dx + 0.5 );
  xnew(wr, n);
  for ( i = 0; i < n; i++ ) {
    x = xmin + (i + 0.5) * dx;
    /* no need to normalize the weight here,
     * it will be done by the WL module */
    wr[i] = pow(x, -exponent);
  }
  return wr;
}



/* initialize HMC state variables,
 * total number of atoms, position, velocity, force */
static void gmxgo_inithmc(gmxgo_t *go, gmx_mtop_t *mtop)
{
  int ib, n = 0;
  gmx_molblock_t *mb = mtop->molblock;
  t_commrec *cr = go->cr;

  /* count the total number of atoms
   * including those the water */
  for ( ib = 0; ib < mtop->nmolblock; ib++ ) {
    n += mb[ib].nmol * mb[ib].natoms_mol;
  }

  if ( DOMAINDECOMP(cr) ) {
    int n1 = (n + cr->nnodes) / cr->nnodes;
    /* we allow some margin */
    n1 = (int) (n1 + 10.0 * sqrt(n1));
    if ( n1 < n ) n = n1;
  }

  go->stncap = (n + gmxgo_blksz - 1) / gmxgo_blksz * gmxgo_blksz;
  go->stn = 0;
  xnew(go->stx, go->stncap);
  xnew(go->stv, go->stncap);
  xnew(go->stf, go->stncap);
  go->hmcrej = 0;
  go->hmctot = 0;
}



/* tp: the temperature
 * ctn: continue from the previous simulation */
static gmxgo_t *gmxgo_open(gmx_mtop_t *mtop, t_commrec *cr,
    double tp, const char *fncfg, int ctn)
{
  gmxgo_t *go;

  xnew(go, 1);

  go->cr = cr;

  /* load input parameters from the configuration file */
  if ( MASTER(cr) ) {
    gmxgocfg_t *cfg = go->cfg;

    gmxgocfg_default(cfg);
    gmxgocfg_load(cfg, fncfg);
  }

  /* broadcast the parameters */
  if ( PAR(cr) ) {
    gmx_bcast(sizeof(gmxgocfg_t), go->cfg, cr);
  }

  go->kT = BOLTZ * tp;


  gmxgo_build(go, mtop, cr);
  if ( PAR(cr) ) {
    go->gvc = gmxvcomm_open(go->n, go->index, cr);
  } else {
    go->gvc = NULL;
  }

  xnew(go->x, go->n);
  xnew(go->f, go->n);
  xnew(go->xref, go->n);
  xnew(go->xf, go->n);
  xnew(go->x1, go->n);
  xnew(go->x2, go->n);
  xnew(go->xwhole, go->n);
  xnew(go->xwholep, go->n);
  xnew(go->xrt, go->n);
  xnew(go->xrtp, go->n);

  /* load the reference */
  if ( MASTER(cr) ) {
    gmxgocfg_t *cfg = go->cfg;

    if ( gmxgo_loadxref(go, cfg->fnpdb) != 0 ) {
      exit(1);
    }

    /* initialize the target distribution */
    go->wr = gmxgo_mkwr(cfg->rmsdmin, cfg->rmsdmax, cfg->rmsddel,
        cfg->wr_exponent);

    /* open a Wang-Landau object */
    go->wl = wl_openf(cfg->rmsdmin, cfg->rmsdmax, cfg->rmsddel,
        cfg->wl_lnf0, cfg->wl_flatness, cfg->wl_frac, cfg->invt_c,
        go->wr, 0);
    if ( go->wl == NULL ) {
      exit(1);
    }

    go->rhis = hist_open(0, cfg->rhis_max, cfg->rhis_dx);

    /* continue from the last run */
    go->isctn = ctn;
    if ( ctn ) {
      wl_load(go->wl, cfg->fnvrmsd);
      hist_load(go->rhis, cfg->fnrhis);
    }

    mtload(cfg->fnmtseed, time(NULL));

    fprintf(stderr, "parallel %d, domain-decomposition %d\n",
        PAR(cr), DOMAINDECOMP(cr));
  }

  gmxgo_inithmc(go, mtop);

  return go;
}



static void gmxgo_close(gmxgo_t *go)
{
  if ( MASTER(go->cr) ) {
    free(go->wr);
    wl_close(go->wl);
    hist_close(go->rhis);
  }

  free(go->index);
  free(go->mass);
  if ( go->gvc != NULL ) {
    gmxvcomm_close(go->gvc);
  }
  free(go->x);
  free(go->f);
  free(go->xref);
  free(go->xf);
  free(go->x1);
  free(go->x2);
  free(go->xwhole);
  free(go->xwholep);
  free(go->xrt);
  free(go->xrtp);
  free(go->stx);
  free(go->stv);
  free(go->stf);
  free(go);
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



/* add the bias force `go->f` to the actual force `f` */
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
/* maximal permissible deviation */
#define CACABOND_DEV 0.06



/* dump basic information
 * including the current RMSD */
__inline static void gmxgo_dump(gmxgo_t *go,
    double (*x)[D], int ePBC, matrix box, gmx_int64_t step)
{
  int i, c1, c;
  double dx[D], dr = 0, f = 0, dev = 0, rmsd;
  t_pbc pbc[1];

  set_pbc(pbc, ePBC, box);

  rmsd = vrmsd(go->xref, go->xf, go->x, go->mass, go->n,
      0, NULL, NULL);
  /* print out the current RMSD */
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



/* save the current position to file
 * `fn0` is the prefix of the file name */
__inline static int gmxgo_writepdb(gmxgo_t *go,
    double (*x)[D], const char *fn0, gmx_int64_t step)
{
  int i;
  char fn[128];
  FILE *fp;

  /* construct a file name with the step number */
  sprintf(fn, "%s_step%s.pdb", fn0, gmx_step_str(step, go->sbuf));

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  for ( i = 0; i < go->n; i++ ) {
    fprintf(fp, "ATOM  %5d  %-3s %3s  %4d    %8.3f%8.3f%8.3f  1.00  1.00          %2c  \n",
        i + 1, go->atnm[i], go->resnm[i], go->resid[i],
        x[i][0] * 10, x[i][1] * 10, x[i][2] * 10, go->atnm[i][0]);
  }

  fclose(fp);
  return 0;
}



/* check the connectivity, that is if all
 * atoms from the same and nearby residues
 * are less than dismax */
static int gmxgo_checkconn(gmxgo_t *go, double (*x)[D],
    double dismax)
{
  int i, j, ir, jr, n = go->n;
  double r2, dis2 = dismax * dismax;

  for ( i = 0; i < n - 1; i++ ) {
    ir = go->resid[i];
    for ( j = i + 1; j < n; j++ ) {
      jr = go->resid[j];
      if ( jr > ir + 1 ) break;
      r2 = vdist2(x[i], x[j]);
      if ( r2 > dis2 ) return -1;
    }
  }
  return 0;
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
    /* we assume that the atoms i and i - 1 are close,
     * and their distance is much smaller than
     * half of the box size.
     * This is true for alpha-carbon atoms.
     * Is it true in general? */
    pbc_dx_d(pbc, xin[i], xin[i - 1], dx);
    if ( go->cfg->rmsdgrp == RMSDGRP_CA ) {
      dr = vnorm(dx);
      dev = fabs(dr - CACABOND);
      if ( dev > devmax ) {
        devmax = dev;
        id = i;
      }
    }
    /* attach the coordinates of i to those of i - 1 */
    vadd(xout[i], xout[i - 1], dx);
  }

  if ( go->cfg->lucky ) return err;

  /* validity tests */
  if ( go->cfg->rmsdgrp == RMSDGRP_CA ) {
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
      hist_save(go->rhis, go->cfg->fnrhis);
      wl_save(go->wl, go->cfg->fnvrmsd);
    }
    //fprintf(stderr, "step %d, max dev %gA\n", (int) step, devmax*10);
  } else {
    /* atoms from the same or nearby residue should not be over 9A apart */
    err = gmxgo_checkconn(go, xout, 0.9);
    if ( err ) {
      fprintf(stderr, "system broken at step %s\n",
          gmx_step_str(step, go->sbuf));
      gmxgo_writepdb(go, xout, "xbadconn", step);
    }
  }

  return err;
}



/* compute the projection of the force `f` along rmsd
 * used in mean-force based flat-histogram sampling
 *
 * NOTE: the mean-force based method is not working
 * so this function is unused now
 * */
static double gmxgo_projectf(gmxgo_t *go, double (*f)[D],
    double rmsd, double kT, double (*x)[D], double (*xf)[D])
{
  int i, n = go->n;
  double dx[D], dot = 0;

  /* the mean force is equal to Sum_i u_i f_i
   * with u_i the conjugate field to the gradient,
   *    dR/dx_i = m_i (x_i - xref_i) / (M R)
   * and u_i should satisfy Sum_i u_i dR/dx_i = 1
   * Here we use u_i = (x_i - xref_i) / R */
  for ( i = 0; i < n; i++ ) {
    /* compute the gradient */
    vdiff(dx, x[i], xf[i]);
    dot += vdot(dx, f[i]);
  }

  /* because of the fixed length of the peptide bonds
   * the factor D * n may not be correct */
  return (dot / kT + D * n - 1) / rmsd;
}



/* compute the force `go->f` of the bias potential
 * for the current `rmsd`
 * compute dV(R)/dR and dV(R)/dx_i, with R being the RMSD
 * and V(R) being the bias potential
 * does not apply to the actual force */
static int gmxgo_calcbiasf(gmxgo_t *go, double rmsd)
{
  int i;
  double dvdx;
  gmxgocfg_t *cfg = go->cfg;

  /* 1. compute dV(R)/dR */
  /* Note that the potential V(x) from the WL module is
   * dimensionless, so we have to multiply it by kT later */
  if ( cfg->bias_mf ) {
    /* compute the gradient from the bias potential */
    dvdx = wl_getdvdx_mf(go->wl, rmsd, cfg->minh,
        cfg->mflmin, cfg->mflmax, cfg->mfhmin, cfg->mfhmax);
  } else {
    /* compute the gradient from the bias mean force
     * Note, the mean-force method is currently unusable
     * only for testing. */
    dvdx = wl_getdvdx_v(go->wl, rmsd,
        cfg->mflmin, cfg->mflmax, cfg->mfhmin, cfg->mfhmax);
  }
  //fprintf(stderr, "rmsd %g, dvdx %g\n", rmsd, dvdx); getchar();

#if 0
  /* the following are debugging code
   * to override the mean force */
  if ( abs(go->cfg->debug) == 15 ) {
    /* In this test, we apply a quadratic potential
     * this test allows us to see if the unit for kT is proper
     * If beta V(R) = 1/2 A (R - R0)^2,
     * then dV(R)/dR = kT A (R - R0),
     * the factor kT is to be applied later on,
     * for a large A, we expect | R - R0 | ~ 1/sqrt(A)
     *
     * Thus, we can compare the distribution with
     *   (0.5*A/pi)**0.5*exp(-0.5*A*(x - R0)**2)
     * in gnuplot
     * */
    double A = 900;
    dvdx = A * (rmsd - 0.5);
  } else if ( abs(go->cfg->debug) == 17 ) {
    /* This potential should roughly cancel
     * this distribution at 600 K
     * which means that the distribution
     * between rmsd = (0.6, 1.0) should be flat */
    if ( fabs(rmsd - 0.8) < 0.2 ) {
      dvdx = -60.0 * (rmsd - 0.8);
    } else {
      /* restraining force */
      dvdx = 100.0 * (rmsd - 0.8);
    }
  }
#endif

  /* multiply the mean force by kB T, where kB is the Boltzmann constant */
  dvdx *= go->kT;

  /* apply the limits */
  if ( dvdx < cfg->mfmin ) {
    dvdx = cfg->mfmin;
  } else if ( dvdx > cfg->mfmax ) {
    dvdx = cfg->mfmax;
  }

  go->dvdx = dvdx;
  //fprintf(stderr, "rmsd %gA, dvdx %g, kT %g (%g, %g)\n", rmsd*10, go->dvdx, go->kT, cfg->mfmin, cfg->mfmax);

  /* 2. compute the force scaling factor
   * f_i = -dV(R)/dx_i
   *     = -(dV(R)/dR) (dR/dx_i)
   *     = -(dV(R)/dR) m_i (x_i - xref_i) / (M R)
   * we divide the factor (M R) first
   **/
  dvdx /= rmsd * go->masstot;

  /* for explicit HMC case, we allow this force
   * only for out-of-boundary cases */
  if ( cfg->exhmc
      && rmsd > go->wl->xmin && rmsd < go->wl->xmax ) {
    dvdx = 0;
    go->dvdx = 0;
  }

  for ( i = 0; i < go->n; i++ ) {
    /* f_i = -[ (dV/dR) / (M R) ] m_i (x_i xref_i) */
    vdiff(go->f[i], go->x[i], go->xf[i]);
    vsmul(go->f[i], -dvdx * go->mass[i]);
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
    go->fplog = fopen(go->cfg->fnlog, (go->isctn ? "a" : "w"));
  }

  if ( go->fplog != NULL ) {
    fprintf(go->fplog, "%s %g\n",
        gmx_step_str(step, go->sbuf), rmsd);
  }
}



/* compute the force from the RMSD bias potential
 * and add this bias force to the actual `f`
 *
 * This function does the following.
 * 1. collect coordinates of the special atom and remove PBC
 * 2. compute the rmsd, and update the bias potential
 *    using the Wang-Landau scheme
 * 3. compute the RMSD bias force
 * 4. add the bias force to the actual force `f`
 *
 * set the argument ePBC as fr->ePBC
 *
 * This routine is called immediate after `do_force()`
 * in `md.c` */
static int gmxgo_rmsd_force(gmxgo_t *go, t_state *state, int doid, rvec *f,
    int ePBC, matrix box, gmx_int64_t step)
{
  double rmsd = 0, mf;
  rvec *x = state->x;
  int err = 0;
  gmxgocfg_t *cfg = go->cfg;
  t_commrec *cr = go->cr;

  /* collect `x` from different nodes to `go->x1` on the master */
  gmxgo_gatherx(go, x, go->x1, doid);

  if ( MASTER(cr) ) {
    /* make the coordinates whole */
    err = gmxgo_shift(go, go->x1, go->x, ePBC, box, step);
    if ( err != 0 ) goto BCAST_ERR;

    //for ( i = 0; i < go->n; i++ )
    //  printf("%d: %g %g %g | %g %g %g\n", i, go->x[i][0], go->x[i][1], go->x[i][2], go->xref[i][0], go->xref[i][1], go->xref[i][2]);

    /* rotate and translate `go->xref` to fit the current
     * configuration `go->x`, the rotated and translated structure
     * is saved in `go->xf`
     * also compute the RMSD between `go->xf` and `go->x` */
    rmsd = vrmsd(go->xref, go->xf, go->x, go->mass, go->n,
        0, NULL, NULL);
    /* save this value, to be reused for explicit HMC in hmcselect */
    go->rmsd = rmsd;

    /* add the current RMSD to WL, update the bias potential */
    wl_addf(go->wl, rmsd);
    /* update the updating magnitude, lnf */
    wl_updatelnf(go->wl);

    /* update the RMSD histogram
     * This histogram is independent of one in WL
     * and it covers a wider range with a finer width */
    hist_add(go->rhis, rmsd);

    if ( (step > 0 && step % cfg->nstchat == 0) || cfg->debug > 0 ) {
      gmxgo_chat(go, step, rmsd);

#if 0
      /* debugging code to test HMC rejection */
      if ( cfg->debug == 3 ) {
        int id;

        for ( id = 0; id < go->stn; id+= 1000 ) {
          fprintf(stderr, "step %d, id %d, x: %g, %g, %g\n", (int) step, id, x[id][0], x[id][1], x[id][2]);
          fprintf(stderr, "step %d, id %d, v: %g, %g, %g\n", (int) step, id, state->v[id][0], state->v[id][1], state->v[id][2]);
          fprintf(stderr, "step %d, id %d, f: %g, %g, %g\n", (int) step, id, f[id][0], f[id][1], f[id][2]);
        }
        //getchar();
      }
#endif

      if ( step > 0 && step % cfg->nstlog == 0 ) {
        gmxgo_log(go, step, rmsd);
      }
    }

    if ( step > 0 && step % cfg->nstrep == 0 ) {
      wl_save(go->wl, cfg->fnvrmsd);
      hist_save(go->rhis, cfg->fnrhis);
      mtsave(cfg->fnmtseed);
    }
  }

BCAST_ERR:
  if ( PAR(cr) ) {
    gmx_barrier(cr);
    gmx_bcast(sizeof(err), &err, cr);
  }
  if ( err ) return err;

  /* compute the projection of force on RMSD
   * only used for the mean-force based histogram flattenning
   * Note, this code is not working, just for testing */
  if ( cfg->bias_mf )
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
    /* compute the force from the RMSD bias potential
     * this call is needed even for explicit HMC in order
     * to handle out-of-boundary cases */
    gmxgo_calcbiasf(go, rmsd);
  }

  /* add the RMSD bias force to the actual force
   * currently this call is made even in the case
   * of explicit HMC in order to handle out of boundary cases
   * the bias force is zero in normal cases, however */
  if ( !cfg->passive ) {
    gmxgo_scatterf(go, f, 0);
  }

  return err;
}



#include "gmx_omp_nthreads.h"

/* copy n vectors */
static void gmxgo_vcopy(rvec *dest, rvec *src, int n)
{
  int th, nth, start_th, end_th;

  nth = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(nth) schedule(static)

  for ( th = 0; th < nth; th++ ) {
    int i, start_th, end_th;
    real *pdest, *psrc;

    start_th = n * DIM * th / nth;
    end_th   = n * DIM * (th + 1) / nth;

    pdest = (real *) dest;
    psrc = (real *) src;

    for ( i = start_th; i < end_th; i++ ) {
      pdest[i] = psrc[i];
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

  if ( go->cfg->passive ) return 0;

  /* determine the number of home atoms */
  n = PAR(cr) ? cr->dd->nat_home : state->natoms;

  if ( n > go->stncap ) {
    int n0 = go->stncap;
    go->stncap = (n + gmxgo_blksz - 1) / gmxgo_blksz * gmxgo_blksz;
    fprintf(stderr, "hmcpushxf: step %s, node %d/%d, expanding HMC state variables %d -> %d (%d)\n",
        gmx_step_str(step, go->sbuf), cr->nodeid, cr->nnodes, n0, go->stncap, n);
    xrenew(go->stx, go->stncap);
    xrenew(go->stv, go->stncap);
    xrenew(go->stf, go->stncap);
  }

#if 0
  if ( go->cfg->debug == 3 && step % 300 - 1 == 3 ) {
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

  if ( go->cfg->passive ) return 0;

  /* determine the number of home atoms */
  n = PAR(cr) ? cr->dd->nat_home : state->natoms;

  /* push v happens after push x, f, so the number of atoms
   * should have been determined */
  if ( go->stn != n ) {
    fprintf(stderr, "step %s, number of atoms mismatch, stn %d, n %d, parallel %d\n",
        gmx_step_str(step, go->sbuf), go->stn, n, PAR(cr));
  }

  gmxgo_vcopy(go->stv, state->v, n);
  return 0;
}



/* pop x, v, f (reverse to the previous configuration)
 * assume the number of home atoms is not changed */
static int gmxgo_hmcpop(gmxgo_t *go,
    t_state *state, rvec *f, int reversev)
{
  int i, n = go->stn, nn, d;
  t_commrec *cr = go->cr;

  if ( go->cfg->passive ) return 0;

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

  return 0;
}



/* compute the RMSD between two structures */
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



/* decide whether to accept the current configuration
 * or return to the previous one
 *
 * This function is called after the GROMACS calls
 * `update_coords()` and `update_constraints()`
 * near the end of the MD step, so the coordinates
 * in `state->x` are the new ones, which differ
 * from the previous coordinates, denoted as `xold`,
 * which are used for
 * `do_force()`,
 * `gmxgo_rmsd_force()`,
 * `gmxgo_hmcpushxf()`.
 *
 * If this function decides to reject the current state,
 * the system will go back to the state of `xold`.
 * */
static int gmxgo_hmcselect(gmxgo_t *go,
    t_state *state, rvec *f, int reversev,
    int ePBC, matrix box, gmx_int64_t step)
{
  double rmsd1, rmsd2, rmsd3, rmsd4, delv, dr0whole, dr0rt;
  int err = 0, acc = 1;
  gmxgocfg_t *cfg = go->cfg;
  t_commrec *cr = go->cr;

  if ( cfg->passive ) return 0;

  /* collect `state->x` from different nodes to `go->x1` on the master */
  gmxgo_gatherx(go, state->x, go->x1, 1);

  if ( MASTER(cr) ) {
    static int once;

    /* remove the periodic boundary condition of `x1`
     * to make the coordinates whole, save the latter in `xwhole` */
    err = gmxgo_shift(go, go->x1, go->xwhole, ePBC, box, step);
    if ( err != 0 ) goto BCAST_ERR;

    /* rotate and translate `xref` to `xrt` to best fit `xwhole`
     * and compute the RMSD between `xrt` (new reference) and `xwhole` (new position) */
    rmsd2 = vrmsd(go->xref, go->xrt, go->xwhole, go->mass, go->n,
        0, NULL, NULL);

    if ( cfg->exhmc ) {
      /* in the explicit HMC scheme, we use the old RMSD
       * from the previous configuration `go->x`,
       * that is, RMSD(go->x, go->xf).
       * which is the configuration saved in `gmxgo_rmsd_force()`
       * right after the `do_force()` call in this step */
      rmsd1 = go->rmsd;
      /*
      rmsd1 = vrmsd(go->xref, NULL, go->x, go->mass, go->n,
          0, NULL, NULL);
      fprintf(stderr, "rmsd %g %g\n", rmsd1, go->rmsd); getchar();
      */
    } else {
      /* compute the RMSD between `xf` (old reference) and `xwhole` (new position)
       * Note that `xf` is essentially the same as `xrtp` except PBC.
       * The periodic cell of the first atom of `xf` should match
       * that of `x` since it is computed in `gmxgo_rmsd_force()`
       * [after `do_force()`, before `gmxgo_hmcpushxf()`]
       * of this step */
      rmsd1 = gmxgo_raw_rmsd(go, go->xf, go->xwhole);

      /* compute the RMSD between `xf` (old reference) and `xwholep` (old position)
       * symmetric to rmsd2 */
      //rmsd4 = go->rmsd;
      rmsd4 = vrmsd(go->xref, go->x1, go->xwholep, go->mass, go->n,
          0, NULL, NULL);

      /* compute the RMSD between `xrt` (new reference) and `xwholep` (old position) */
      rmsd3 = gmxgo_raw_rmsd(go, go->xrt, go->xwholep);
    }

    /* checking code to see if our assumption on the PBC is correct */
    if ( !cfg->lucky && once ) {
      /* the two whole structures `x` and `xwhole` are supposed to be close,
       * if not, one of them might have been wrapped unintendedly
       * Note `xwholep` is supposed to be same as `x`, the time evolution is
       *   xwholep --> x --> x1 --> xwhole
       * where in the first step, x might have been wrapped because of
       * the periodic boundary condition */
      dr0whole = vdist(go->x[0], go->xwhole[0]);
      /* the two fitted structures `xf` and `xrt` are supposed to be close
       *   xrpt --> xf --> xrt */
      if ( cfg->exhmc ) {
        dr0rt = 0;
      } else {
        dr0rt = vdist(go->xf[0], go->xrt[0]);
      }

      if ( dr0whole > 0.1 || dr0rt > 0.5 ) {
        /* the threshold 0.1nm or 1A is generous, can be smaller */
        fprintf(stderr, "\nWarning. step %s: the first special atoms of the two "
            "whole/fitted structures differ by %gA (whole), %gA (fit)\n",
            gmx_step_str(step, go->sbuf), dr0whole * 10, dr0rt * 10);
        /* xwholep == x, xrtp == xf */
        gmxgo_writepdb(go, go->xwholep, "xwholep", step);
        gmxgo_writepdb(go, go->xwhole,  "xwhole",  step);
        gmxgo_writepdb(go, go->xrtp,    "xrtp",    step);
        gmxgo_writepdb(go, go->xrt,     "xrt",     step);
        gmxgo_writepdb(go, go->xf,      "xf",      step);
        gmxgo_writepdb(go, go->x,       "x",       step);
        gmxgo_writepdb(go, go->x1,      "x1",      step);
      }
    }

    if ( cfg->bias_mf ) {
      /* this branch is not currently used */
      delv = wl_getdelv_mf(go->wl, rmsd1, rmsd2, cfg->minh,
          cfg->mflmin, cfg->mflmax, cfg->mfhmin, cfg->mfhmax);
    } else {
      delv = wl_getdelv_v(go->wl, rmsd1, rmsd2);
      if ( !cfg->exhmc ) {
        /* for implicit HMC */
        if ( once ) {
          delv -= wl_getdelv_v(go->wl, rmsd3, rmsd4);
          delv *= 0.5;
        }

        /* we don't want trouble for boundary cases */
        if ( ( rmsd1 <= go->wl->xmin || rmsd1 >= go->wl->xmax )
          || ( rmsd2 <= go->wl->xmin || rmsd2 >= go->wl->xmax )
          || ( rmsd3 <= go->wl->xmin || rmsd3 >= go->wl->xmax )
          || ( rmsd4 <= go->wl->xmin || rmsd4 >= go->wl->xmax ) ) {
          delv = 0;
        }
      }
    }

    acc = 1;
    if ( delv > 0 ) {
      double r = rand01();
      acc = ( r < exp(-delv) );
    }

#if 0
    /* special flag to test HMC rejection */
    if ( cfg->debug == 3 ) {
      acc = ( step % 300 != 3 );
    }
#endif

    if ( !acc ) {
      go->hmcrej += 1;
      if ( !cfg->exhmc ) {
        /* it is rare that implicit HMC rejects a state */
        fprintf(stderr, "step %s: HMC rejecting a state, rmsd %g -> %g, [3: %g -> 4: %g] delv %g\n",
            gmx_step_str(step, go->sbuf), rmsd1, rmsd2, rmsd3, rmsd4, delv);
      }
    }
    go->hmctot += 1;

    /* backup the structures */
    memcpy(go->xwholep, go->xwhole, go->n * sizeof(go->xwhole[0]));
    if ( !cfg->exhmc ) {
      memcpy(go->xrtp, go->xrt, go->n * sizeof(go->xrt[0]));
    }
    once = 1;

    /* TODO: we may need to scramble velocities for explicit HMC */
  }

BCAST_ERR:
  if ( PAR(cr) ) {
    gmx_bcast(sizeof(err), &err, cr);
    gmx_bcast(sizeof(acc), &acc, cr);
  }

  if ( !err && !acc ) {
    gmxgo_hmcpop(go, state, f, reversev);
  }

  if ( !cfg->lucky && PAR(cr) ) {
    gmx_barrier(cr);
  }

  return acc;
}



#endif /* GMXGO_H__ */
