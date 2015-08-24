#ifndef HMCRMSD_H__
#define HMCRMSD_H__

#include "util.h"



typedef struct {
  const char *fnpdb;
} hmcrmsd_t;



/* load settings from the configuration file `fn' */
static int hmcrmsd_loadcfg(hmcrmsd_t *h, const char *fn)
{
  FILE *fp;
  char buf[800], *p, *key, *val;
  int inpar;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return -1;
  }

  while ( fgets(buf, sizeof buf, fp) ) {
    strstrip(buf); /* remove trailing spaces */
    if ( buf[0] == '\0' || buf[0] == '#' ) continue;

    /* parse the line to a key and a value */
    /* find the end of the key */
    inpar = 0; /* within (...) */
    for ( p = buf; *p; p++ ) {
      if ( (!inpar && isspace(*p)) || *p == '=' ) {
        *p = '\0'; /* end the key part */
        break;
      }
      if ( !inpar && (*p == '(' || *p == '[') )
        inpar = 1; /* enter a parentheses block */
      else if ( inpar && (*p == ')' || *p == ']') )
        inpar = 0; /* leave a parentheses block */
      *p = (char) tolower(*p);
    }
    key = buf;

    /* find the beginning of the value */
    for ( p++; isspace(*p) || *p == '=' ; ) p++;
    val = p;

    if ( strcmpfuzzy(key, "pdb") == 0
      || strcmpfuzzy(key, "fnpdb") == 0 ) {
      h->fnpdb = val;
    } else {
      fprintf(stderr, "Warning: unknown option %s = %s\n",
          key, val);
      getchar();
    }
  }

  fclose(fp);

  return 0;
}



/* open an HMC object for a flat histogram sampling
 * along the RMSD */
static hmcrmsd_t *hmcrmsd_open(const char *fncfg)
{
  hmcrmsd_t *h;

  if ((h = (hmcrmsd_t *) calloc(1, sizeof(*h))) == NULL ) {
    fprintf(stderr, "no memory for HMCRMSD\n");
    return NULL;
  }

  /* set the default values */
  h->fnpdb = "ref.pdb";

  fprintf(stderr, "loading configuration from %s\n", fncfg);

  if ( hmcrmsd_loadcfg(h, fncfg) != 0 ) {
    exit(1);
    return NULL;
  }

  fprintf(stderr, "PDB %s\n", h->fnpdb);

  return h;
}



static void hmcrmsd_close(hmcrmsd_t *h)
{
  free(h);
}



#endif /* HMCRMSD_H__ */
