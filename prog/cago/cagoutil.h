#ifndef CAGOUTIL_H__
#define CAGOUTIL_H__



#include "mdutil.h"
#include "mtrand.h"
#include "util.h"



const char *aanames[] = {
  "GLY", "ALA", "VAL", "LEU", "ILE",
  "PRO", "SER", "THR", "CYS", "MET",
  "ASN", "GLN", "ASP", "GLU",
  "LYS", "ARG", "HIS",
  "PHE", "TYR", "TRP"};

const char aaletters[] = "GAVLIPSTCMNQDEKRHFYW";

/* residue name to integer of amino acid */
__inline static int res2iaa(const char *res)
{
  int i;

  for ( i = 0; i < 20; i++ ) {
    if ( strcmp(res, aanames[i]) == 0 ) {
      return i;
    }
  }
  return -1;
}



#endif /* CAGOUTIL_H__ */

