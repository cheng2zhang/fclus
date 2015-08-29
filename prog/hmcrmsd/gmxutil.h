#ifndef GMXUTIL_H__
#define GMXUTIL_H__


/* Debugging utilities */

/* copy the definition of t_gmx_update here to access upd->xp
 * which is defined in gromacs/mdlib/update.c */
typedef struct gmx_stochd_t *gmx_stochd_t_ptr;

typedef struct gmx_update {
    gmx_stochd_t_ptr sd;
    rvec         *xp;
    int           xp_nalloc;
    gmx_int64_t     deformref_step;
    matrix          deformref_box;
} t_gmx_update;



#endif /* defined(GMXUTIL_H__) */
