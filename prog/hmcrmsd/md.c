/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "sysstuff.h"
#include "vec.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "calcmu.h"
#include "index.h"
#include "update.h"
#include "ns.h"
#include "mdrun.h"
#include "md_support.h"
#include "md_logging.h"
#include "network.h"
#include "xvgr.h"
#include "physics.h"
#include "names.h"
#include "force.h"
#include "disre.h"
#include "orires.h"
#include "pme.h"
#include "mdatoms.h"
#include "deform.h"
#include "qmmm.h"
#include "domdec.h"
#include "domdec_network.h"
#include "gromacs/gmxlib/topsort.h"
#include "coulomb.h"
#include "constr.h"
#include "gromacs/gmxpreprocess/compute_io.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "sighandler.h"
#include "txtdump.h"
#include "gromacs/utility/cstringutil.h"
#include "pme_loadbal.h"
#include "bondf.h"
#include "types/nlistheuristics.h"
#include "types/iteratedconstraints.h"
#include "nbnxn_cuda_data_mgmt.h"

#include "gromacs/utility/gmxmpi.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trajectory_writing.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/imd/imd.h"

#include "gmxgo.h"

static void reset_all_counters(FILE *fplog, t_commrec *cr,
                               gmx_int64_t step,
                               gmx_int64_t *step_rel, t_inputrec *ir,
                               gmx_wallcycle_t wcycle, t_nrnb *nrnb,
                               gmx_walltime_accounting_t walltime_accounting,
                               nbnxn_cuda_ptr_t cu_nbv)
{
    char sbuf[STEPSTRSIZE];

    /* Reset all the counters related to performance over the run */
    md_print_warn(cr, fplog, "step %s: resetting all time and cycle counters\n",
                  gmx_step_str(step, sbuf));

    if (cu_nbv)
    {
        nbnxn_cuda_reset_timings(cu_nbv);
    }

    wallcycle_stop(wcycle, ewcRUN);
    wallcycle_reset_all(wcycle);
    if (DOMAINDECOMP(cr))
    {
        reset_dd_statistics_counters(cr->dd);
    }
    init_nrnb(nrnb);
    ir->init_step += *step_rel;
    ir->nsteps    -= *step_rel;
    *step_rel      = 0;
    wallcycle_start(wcycle, ewcRUN);
    walltime_accounting_start(walltime_accounting);
    print_date_and_time(fplog, cr->nodeid, "Restarted time", gmx_gettime());
}

double mdhmcrmsd(FILE *fplog, t_commrec *cr, int nfile, const t_filenm fnm[],
             const output_env_t oenv, gmx_bool bVerbose, gmx_bool bCompact,
             int nstglobalcomm,
             gmx_constr_t constr,
             int stepout, t_inputrec *ir,
             gmx_mtop_t *top_global,
             t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb, gmx_wallcycle_t wcycle,
             gmx_edsam_t ed, t_forcerec *fr,
             real cpt_period, real max_hours,
             const char gmx_unused *deviceOptions,
             int imdport,
             unsigned long Flags,
             gmx_walltime_accounting_t walltime_accounting)
{
    gmx_mdoutf_t    outf = NULL;
    gmx_int64_t     step, step_rel;
    double          elapsed_time;
    double          t, t0, lam0[efptNR];
    gmx_bool        bGStatEveryStep, bGStat, bCalcVir, bCalcEner;
    gmx_bool        bNS, bNStList, bSimAnn, bStopCM, bNotLastFrame = FALSE,
                    bFirstStep, bStateFromCP, bStateFromTPX, bInitStep, bLastStep,
                    bBornRadii, bStartingFromCpt;
    gmx_bool          do_ene, do_log, do_verbose,
                      bForceUpdate = FALSE, bCPT;
    gmx_bool          bMasterState;
    int               force_flags, cglo_flags;
    tensor            force_vir, shake_vir, total_vir, tmp_vir, pres;
    int               i, m;
    t_trxstatus      *status;
    rvec              mu_tot;
    t_vcm            *vcm;
    t_state          *bufstate = NULL;
    matrix           *scale_tot, pcoupl_mu, M, ebox;
    gmx_nlheur_t      nlh;
    int               nchkpt  = 1;
    gmx_localtop_t   *top;
    t_mdebin         *mdebin   = NULL;
    t_state          *state    = NULL;
    rvec             *f_global = NULL;
    gmx_enerdata_t   *enerd;
    rvec             *f = NULL;
    gmx_global_stat_t gstat;
    gmx_update_t      upd   = NULL;
    t_graph          *graph = NULL;
    globsig_t         gs;
    gmx_groups_t     *groups;
    gmx_ekindata_t   *ekind, *ekind_save;
    int               count, nconverged = 0;
    real              timestep   = 0;
    double            tcount     = 0;
    gmx_bool          bConverged = TRUE, bOK, bSumEkinhOld, bNeedRepartition;
    gmx_bool          bAppend;
    gmx_bool          bResetCountersHalfMaxH = FALSE;
    gmx_bool          bVV, bIterativeCase, bFirstIterate, bTemp, bPres, bTrotter;
    gmx_bool          bUpdateDoLR;
    real              dvdl_constr;
    rvec             *cbuf        = NULL;
    int               cbuf_nalloc = 0;
    matrix            lastbox;
    real              veta_save, scalevir, tracevir;
    real              vetanew = 0;
    int               lamnew  = 0;
    /* for FEP */
    double            cycles;
    real              saved_conserved_quantity = 0;
    real              last_ekin                = 0;
    int               iter_i;
    t_extmass         MassQ;
    int             **trotter_seq;
    char              sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
    int               handled_stop_condition = gmx_stop_cond_none; /* compare to get_stop_condition*/
    gmx_iterate_t     iterate;
    gmx_int64_t       multisim_nsteps = -1;                        /* number of steps to do  before first multisim
                                                                          simulation stops. If equal to zero, don't
                                                                          communicate any more between multisims.*/
    /* PME load balancing data for GPU kernels */
    pme_load_balancing_t pme_loadbal = NULL;
    double               cycles_pmes;
    gmx_bool             bPMETuneTry = FALSE, bPMETuneRunning = FALSE;

    /* Interactive MD */
    gmx_bool          bIMDstep = FALSE;

    /* object of handling the bias potential based on RMSD */
    gmxgo_t *go;


    /* Check for special mdrun options */
    bAppend  = (Flags & MD_APPENDFILES);
    if (Flags & MD_RESETCOUNTERSHALFWAY)
    {
        if (ir->nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle, ir->nsteps/2);
        }
        /* Signal to reset the counters halfway the simulation time. */
        bResetCountersHalfMaxH = (max_hours > 0);
    }

    /* md-vv uses averaged full step velocities for T-control
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bVV = EI_VV(ir->eI);
    /* all the iteratative cases - only if there are constraints */
    bIterativeCase = ((IR_NPH_TROTTER(ir) || IR_NPT_TROTTER(ir)) && (constr));
    gmx_iterate_init(&iterate, FALSE); /* The default value of iterate->bIterationActive is set to
                                          false in this step.  The correct value, true or false,
                                          is set at each step, as it depends on the frequency of temperature
                                          and pressure control.*/
    bTrotter = (bVV && (IR_NPT_TROTTER(ir) || IR_NPH_TROTTER(ir) || IR_NVT_TROTTER(ir)));

    check_ir_old_tpx_versions(cr, fplog, ir, top_global);

    nstglobalcomm   = check_nstglobalcomm(fplog, cr, nstglobalcomm, ir);
    bGStatEveryStep = (nstglobalcomm == 1);

    if (!bGStatEveryStep && ir->nstlist == -1 && fplog != NULL)
    {
        fprintf(fplog,
                "To reduce the energy communication with nstlist = -1\n"
                "the neighbor list validity should not be checked at every step,\n"
                "this means that exact integration is not guaranteed.\n"
                "The neighbor list validity is checked after:\n"
                "  <n.list life time> - 2*std.dev.(n.list life time)  steps.\n"
                "In most cases this will result in exact integration.\n"
                "This reduces the energy communication by a factor of 2 to 3.\n"
                "If you want less energy communication, set nstlist > 3.\n\n");
    }

    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog, cr, ir, oenv, &t, &t0, state_global->lambda,
            &(state_global->fep_state), lam0,
            nrnb, top_global, &upd,
            nfile, fnm, &outf, &mdebin,
            force_vir, shake_vir, mu_tot, &bSimAnn, &vcm, Flags, wcycle);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd, 1);
    init_enerdata(top_global->groups.grps[egcENER].nr, ir->fepvals->n_lambda,
                  enerd);
    if (DOMAINDECOMP(cr))
    {
        f = NULL;
    }
    else
    {
        snew(f, top_global->natoms);
    }

    /* Kinetic energy data */
    snew(ekind, 1);
    init_ekindata(fplog, top_global, &(ir->opts), ekind);
    /* needed for iteration of constraints */
    snew(ekind_save, 1);
    init_ekindata(fplog, top_global, &(ir->opts), ekind_save);
    /* Copy the cos acceleration to the groups struct */
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);
    debug_gmx();

    if (DEFORM(*ir))
    {
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
        set_deform_reference_box(upd,
                                 deform_init_init_step_tpx,
                                 deform_init_box_tpx);
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
    }

    {
        double io = compute_io(ir, top_global->natoms, groups, mdebin->ebin->nener, 1);
        if ((io > 2000) && MASTER(cr))
        {
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        top = dd_init_local_top(top_global);

        snew(state, 1);
        dd_init_local_state(cr->dd, state_global, state);

        if (DDMASTER(cr->dd) && ir->nstfout)
        {
            snew(f_global, state_global->natoms);
        }
    }
    else
    {
        top = gmx_mtop_generate_local_top(top_global, ir);

        forcerec_set_excl_load(fr, top);

        state    = serial_init_local_state(state_global);
        f_global = f;

        atoms2md(top_global, ir, 0, NULL, top_global->natoms, mdatoms);

        if (ir->ePBC != epbcNONE && !fr->bMolPBC)
        {
            graph = mk_graph(fplog, &(top->idef), 0, top_global->natoms, FALSE, FALSE);
        }

        setup_bonded_threading(fr, &top->idef);
    }

    /* Set up interactive MD (IMD) */
    init_IMD(ir, cr, top_global, fplog, ir->nstcalcenergy, state_global->x,
             nfile, fnm, oenv, imdport, Flags);

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, ir->init_step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            state, &f, mdatoms, top, fr,
                            NULL, NULL, constr,
                            nrnb, wcycle, FALSE);

    }

    update_mdatoms(mdatoms, state->lambda[efptMASS]);

    if (opt2bSet("-cpi", nfile, fnm))
    {
        bStateFromCP = gmx_fexist_master(opt2fn_master("-cpi", nfile, fnm, cr), cr);
    }
    else
    {
        bStateFromCP = FALSE;
    }

    /* initialize an object for Go model
     * all nodes call this */
    go = gmxgo_open(top_global, cr, ir->opts.ref_t[0],
        opt2fn_master("-cfg", nfile, fnm, cr), bStateFromCP);

    if (MASTER(cr))
    {
        if (bStateFromCP)
        {
            /* Update mdebin with energy history if appending to output files */
            if (Flags & MD_APPENDFILES)
            {
                restore_energyhistory_from_state(mdebin, &state_global->enerhist);
            }
            else
            {
                /* We might have read an energy history from checkpoint,
                 * free the allocated memory and reset the counts.
                 */
                done_energyhistory(&state_global->enerhist);
                init_energyhistory(&state_global->enerhist);
            }
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(&state_global->enerhist, mdebin);
    }

    /* Initialize constraints */
    if (constr && !DOMAINDECOMP(cr))
    {
        set_constraints(constr, top, ir, mdatoms, cr);
    }

    /* PME tuning is only supported with GPUs or PME nodes.
     * PME tuning is not supported with PME only for LJ and not for Coulomb.
     */
    if ((Flags & MD_TUNEPME) &&
        EEL_PME(fr->eeltype) &&
        ( (fr->cutoff_scheme == ecutsVERLET && fr->nbv->bUseGPU) || !(cr->duty & DUTY_PME)))
    {
        pme_loadbal_init(&pme_loadbal, ir, state->box, fr->ic, fr->pmedata);
        cycles_pmes = 0;
        if (cr->duty & DUTY_PME)
        {
            /* Start tuning right away, as we can't measure the load */
            bPMETuneRunning = TRUE;
        }
        else
        {
            /* Separate PME nodes, we can measure the PP/PME load balance */
            bPMETuneTry = TRUE;
        }
    }

    if (!ir->bContinuation)
    {
        if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
        {
            /* Set the velocities of frozen particles to zero */
            for (i = 0; i < mdatoms->homenr; i++)
            {
                for (m = 0; m < DIM; m++)
                {
                    if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m])
                    {
                        state->v[i][m] = 0;
                    }
                }
            }
        }

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog, constr, ir, mdatoms, state,
                               cr, nrnb, fr, top);
        }
    }

    debug_gmx();

    /* Be REALLY careful about what flags you set here. You CANNOT assume
     * this is the first step, since we might be restarting from a checkpoint,
     * and in that case we should not do any modifications to the state.
     */
    bStopCM = (ir->comm_mode != ecmNO && !ir->bContinuation);

    cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
                  | (bStopCM ? CGLO_STOPCM : 0)
                  | (bVV ? CGLO_PRESSURE : 0)
                  | (bVV ? CGLO_CONSTRAINT : 0)
                  | ((Flags & MD_READ_EKIN) ? CGLO_READEKIN : 0));

    bSumEkinhOld = FALSE;
    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                    NULL, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                    constr, NULL, FALSE, state->box,
                    top_global, &bSumEkinhOld, cglo_flags);
    if (ir->eI == eiVVAK)
    {
        /* a second call to get the half step temperature initialized as well */
        /* we do the same call as above, but turn the pressure off -- internally to
           compute_globals, this is recognized as a velocity verlet half-step
           kinetic energy calculation.  This minimized excess variables, but
           perhaps loses some logic?*/

        compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                        NULL, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                        constr, NULL, FALSE, state->box,
                        top_global, &bSumEkinhOld,
                        cglo_flags &~(CGLO_STOPCM | CGLO_PRESSURE));
    }

    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!(Flags & MD_STARTFROMCPT))
    {
        for (i = 0; (i < ir->opts.ngtc); i++)
        {
            copy_mat(ekind->tcstat[i].ekinh, ekind->tcstat[i].ekinh_old);
        }
    }
    if (ir->eI != eiVV)
    {
        enerd->term[F_TEMP] *= 2; /* result of averages being done over previous and current step,
                                     and there is no previous step */
    }

    /* if using an iterative algorithm, we need to create a working directory for the state. */
    if (bIterativeCase)
    {
        bufstate = init_bufstate(state);
    }

    /* need to make an initiation call to get the Trotter variables set, as well as other constants for non-trotter
       temperature control */
    trotter_seq = init_npt_vars(ir, state, &MassQ, bTrotter);

    if (MASTER(cr))
    {
        if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
        {
            fprintf(fplog,
                    "RMS relative constraint deviation after constraining: %.2e\n",
                    constr_rmsd(constr, FALSE));
        }
        if (EI_STATE_VELOCITY(ir->eI))
        {
            fprintf(fplog, "Initial temperature: %g K\n", enerd->term[F_TEMP]);
        }
        {
            char tbuf[20];
            fprintf(stderr, "starting mdrun '%s'\n",
                    *(top_global->name));
            if (ir->nsteps >= 0)
            {
                sprintf(tbuf, "%8.1f", (ir->init_step+ir->nsteps)*ir->delta_t);
            }
            else
            {
                sprintf(tbuf, "%s", "infinite");
            }
            if (ir->init_step > 0)
            {
                fprintf(stderr, "%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                        gmx_step_str(ir->init_step+ir->nsteps, sbuf), tbuf,
                        gmx_step_str(ir->init_step, sbuf2),
                        ir->init_step*ir->delta_t);
            }
            else
            {
                fprintf(stderr, "%s steps, %s ps.\n",
                        gmx_step_str(ir->nsteps, sbuf), tbuf);
            }
        }
        fprintf(fplog, "\n");
    }

    walltime_accounting_start(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, "mdrun");

    /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
    debug_gmx();
    /***********************************************************
     *
     *             Loop over MD steps
     *
     ************************************************************/

    /* loop over MD steps */
    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bStateFromTPX    = !bStateFromCP;
    bInitStep        = bFirstStep && (bStateFromTPX || bVV);
    bStartingFromCpt = (Flags & MD_STARTFROMCPT) && bInitStep;
    bLastStep        = FALSE;
    bSumEkinhOld     = FALSE;
    bNeedRepartition = FALSE;

    init_global_signals(&gs, cr, ir, 0);

    step     = ir->init_step;
    step_rel = 0;

    if (ir->nstlist == -1)
    {
        init_nlistheuristics(&nlh, bGStatEveryStep, step);
    }

    if (MULTISIM(cr))
    {
        /* check how many steps are left in other sims */
        multisim_nsteps = get_multisim_nsteps(cr, ir->nsteps);
    }

    /* CZ: force temperature coupling at every step
     * this may be unnecessary, but it is required by reversibility
     * when using VV, we probably wish for a reversible trajectory
     * or set it to 1.0 helps */
    if ( EI_VV(ir->eI) ) {
      ir->nsttcouple = 1;
    }

    /* and stop now if we should */
    bLastStep = ((ir->nsteps >= 0 && step_rel > ir->nsteps) ||
                 ((multisim_nsteps >= 0) && (step_rel >= multisim_nsteps )));
    while (!bLastStep)
    {

        wallcycle_start(wcycle, ewcSTEP);

        {
            bLastStep = (step_rel == ir->nsteps);
            t         = t0 + step*ir->delta_t;
        }

        if (bSimAnn)
        {
            update_annealing_target_temp(&(ir->opts), t);
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

        {
            /* Determine whether or not to do Neighbour Searching and LR */
            bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);

            bNS = (bFirstStep || bNeedRepartition || bNStList ||
                   (ir->nstlist == -1 && nlh.nabnsb > 0));

            if (bNS && ir->nstlist == -1)
            {
                set_nlistheuristics(&nlh, bFirstStep || bNeedRepartition, step);
            }
        }

        /* check whether we should stop because another simulation has
           stopped. */
        if (MULTISIM(cr))
        {
            if ( (multisim_nsteps >= 0) &&  (step_rel >= multisim_nsteps)  &&
                 (multisim_nsteps != ir->nsteps) )
            {
                if (bNS)
                {
                    if (MASTER(cr))
                    {
                        fprintf(stderr,
                                "Stopping simulation %d because another one has finished\n",
                                cr->ms->sim);
                    }
                    bLastStep         = TRUE;
                    gs.sig[eglsCHKPT] = 1;
                }
            }
        }

        /* < 0 means stop at next step, > 0 means stop at next NS step */
        if ( (gs.set[eglsSTOPCOND] < 0) ||
             ( (gs.set[eglsSTOPCOND] > 0) && (bNStList || ir->nstlist == 0) ) )
        {
            bLastStep = TRUE;
        }

        /* Determine whether or not to update the Born radii if doing GB */
        bBornRadii = bFirstStep;
        if (ir->implicit_solvent && (step % ir->nstgbradii == 0))
        {
            bBornRadii = TRUE;
        }

        do_log     = do_per_step(step, ir->nstlog) || bFirstStep || bLastStep;
        do_verbose = bVerbose &&
            (step % stepout == 0 || bFirstStep || bLastStep);

        if (bNS && !(bFirstStep && ir->bContinuation))
        {
            {
                bMasterState = FALSE;
                /* Correct the new box if it is too skewed */
                if (DYNAMIC_BOX(*ir))
                {
                    if (correct_box(fplog, step, state->box, graph))
                    {
                        bMasterState = TRUE;
                    }
                }
                if (DOMAINDECOMP(cr) && bMasterState)
                {
                    dd_collect_state(cr->dd, state, state_global);
                }
            }

            if (DOMAINDECOMP(cr))
            {
                /* Repartition the domain decomposition */
                wallcycle_start(wcycle, ewcDOMDEC);
                dd_partition_system(fplog, step, cr,
                                    bMasterState, nstglobalcomm,
                                    state_global, top_global, ir,
                                    state, &f, mdatoms, top, fr,
                                    NULL, NULL, constr,
                                    nrnb, wcycle,
                                    do_verbose && !bPMETuneRunning);
                wallcycle_stop(wcycle, ewcDOMDEC);
                /* If using an iterative integrator, reallocate space to match the decomposition */
            }
        }

        if (MASTER(cr) && do_log)
        {
            print_ebin_header(fplog, step, t, state->lambda[efptFEP]); /* can we improve the information printed here? */
        }

        clear_mat(force_vir);

        /* We write a checkpoint at this MD step when:
         * either at an NS step when we signalled through gs,
         * or at the last step (but not when we do not want confout),
         * but never at the first step.
         */
        bCPT = (((gs.set[eglsCHKPT] && (bNS || ir->nstlist == 0)) ||
                 (bLastStep && (Flags & MD_CONFOUT))) &&
                step > ir->init_step);
        if (bCPT)
        {
            gs.set[eglsCHKPT] = 0;
        }

        /* Determine the energy and pressure:
         * at nstcalcenergy steps and at energy output steps (set below).
         */
        if (EI_VV(ir->eI) && (!bInitStep))
        {
            /* for vv, the first half of the integration actually corresponds
               to the previous step.  bCalcEner is only required to be evaluated on the 'next' step,
               but the virial needs to be calculated on both the current step and the 'next' step. Future
               reorganization may be able to get rid of one of the bCalcVir=TRUE steps. */

            bCalcEner = do_per_step(step-1, ir->nstcalcenergy);
            bCalcVir  = bCalcEner ||
                (ir->epc != epcNO && (do_per_step(step, ir->nstpcouple) || do_per_step(step-1, ir->nstpcouple)));
        }
        else
        {
            bCalcEner = do_per_step(step, ir->nstcalcenergy);
            bCalcVir  = bCalcEner ||
                (ir->epc != epcNO && do_per_step(step, ir->nstpcouple));
        }

        /* Do we need global communication ? */
        bGStat = (bCalcVir || bCalcEner || bStopCM ||
                  do_per_step(step, nstglobalcomm) || (bVV && IR_NVT_TROTTER(ir) && do_per_step(step-1, nstglobalcomm)) ||
                  (ir->nstlist == -1 && step >= nlh.step_nscheck));

        do_ene = (do_per_step(step, ir->nstenergy) || bLastStep);

        if (do_ene || do_log)
        {
            bCalcVir  = TRUE;
            bCalcEner = TRUE;
            bGStat    = TRUE;
        }

        /* these CGLO_ options remain the same throughout the iteration */
        cglo_flags = (
                      (bGStat ? CGLO_GSTAT : 0)
                      );

        force_flags = (GMX_FORCE_STATECHANGED |
                       (DYNAMIC_BOX(*ir) ? GMX_FORCE_DYNAMICBOX : 0) |
                       GMX_FORCE_ALLFORCES |
                       GMX_FORCE_SEPLRF |
                       (bCalcVir ? GMX_FORCE_VIRIAL : 0) |
                       (bCalcEner ? GMX_FORCE_ENERGY : 0)
                       );

        if (fr->bTwinRange)
        {
            if (do_per_step(step, ir->nstcalclr))
            {
                force_flags |= GMX_FORCE_DO_LR;
            }
        }

        {
            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */
            do_force(fplog, cr, ir, step, nrnb, wcycle, top, groups,
                     state->box, state->x, &state->hist,
                     f, force_vir, mdatoms, enerd, fcd,
                     state->lambda, graph,
                     fr, NULL, mu_tot, t, mdoutf_get_fp_field(outf), ed, bBornRadii,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags);
        }

        /* compute the bias force */
        if ( 0 !=  gmxgo_rmsd_force(go, state, 1, f, fr->ePBC, state->box, step) ) {
          fprintf(stderr, "step %s: gmxgo_rmsd_force failed\n", gmx_step_str(step, go->sbuf));
          exit(1);
        }

        /* push coordinates and force for HMC
         * In the VV case, we are still in previous half-step
         * so the pushed coordinates and force are those
         * at the beginning of the current step
         * In the leapfrog case, this is done before the
         * velocity update */
        if ( go->cfg->dohmc ) {
          gmxgo_hmcpushxf(go, state, f, step);
        }

        if (bVV && !bStartingFromCpt)
        /*  ############### START FIRST UPDATE HALF-STEP FOR VV METHODS############### */
        {
            rvec *vbuf = NULL;

            wallcycle_start(wcycle, ewcUPDATE);
            if (ir->eI == eiVV && bInitStep)
            {
                /* if using velocity verlet with full time step Ekin,
                 * take the first half step only to compute the
                 * virial for the first step. From there,
                 * revert back to the initial coordinates
                 * so that the input is actually the initial step.
                 */
                snew(vbuf, state->natoms);
                copy_rvecn(state->v, vbuf, 0, state->natoms); /* should make this better for parallelizing? */
            }
            else
            {
                /* this is for NHC in the Ekin(t+dt/2) version of vv */
                trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ1);
            }

            /* If we are using twin-range interactions where the long-range component
             * is only evaluated every nstcalclr>1 steps, we should do a special update
             * step to combine the long-range forces on these steps.
             * For nstcalclr=1 this is not done, since the forces would have been added
             * directly to the short-range forces already.
             *
             * TODO Remove various aspects of VV+twin-range in master
             * branch, because VV integrators did not ever support
             * twin-range multiple time stepping with constraints.
             */
            bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

            update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC,
                          f, bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                          ekind, M, upd, bInitStep, etrtVELOCITY1,
                          cr, nrnb, constr, &top->idef);

            if (bIterativeCase && do_per_step(step-1, ir->nstpcouple) && !bInitStep)
            {
                gmx_iterate_init(&iterate, TRUE);
            }
            /* for iterations, we save these vectors, as we will be self-consistently iterating
               the calculations */

            /*#### UPDATE EXTENDED VARIABLES IN TROTTER FORMULATION */

            /* save the state */
            if (iterate.bIterationActive)
            {
                copy_coupling_state(state, bufstate, ekind, ekind_save, &(ir->opts));
            }

            bFirstIterate = TRUE;
            while (bFirstIterate || iterate.bIterationActive)
            {
                if (iterate.bIterationActive)
                {
                    copy_coupling_state(bufstate, state, ekind_save, ekind, &(ir->opts));
                    if (bFirstIterate && bTrotter)
                    {
                        /* The first time through, we need a decent first estimate
                           of veta(t+dt) to compute the constraints.  Do
                           this by computing the box volume part of the
                           trotter integration at this time. Nothing else
                           should be changed by this routine here.  If
                           !(first time), we start with the previous value
                           of veta.  */

                        veta_save = state->veta;
                        trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ0);
                        vetanew     = state->veta;
                        state->veta = veta_save;
                    }
                }

                bOK = TRUE;
                {
                    wallcycle_stop(wcycle, ewcUPDATE);
                    update_constraints(fplog, step, NULL, ir, ekind, mdatoms,
                                       state, fr->bMolPBC, graph, f,
                                       &top->idef, shake_vir,
                                       cr, nrnb, wcycle, upd, constr,
                                       TRUE, bCalcVir, vetanew);
                    wallcycle_start(wcycle, ewcUPDATE);

                    if (bCalcVir && bUpdateDoLR && ir->nstcalclr > 1)
                    {
                        /* Correct the virial for multiple time stepping */
                        m_sub(shake_vir, fr->vir_twin_constr, shake_vir);
                    }

                    if (!bOK)
                    {
                        gmx_fatal(FARGS, "Constraint error: Shake, Lincs or Settle could not solve the constrains");
                    }

                }

                /* if VV, compute the pressure and constraints */
                /* For VV2, we strictly only need this if using pressure
                 * control, but we really would like to have accurate pressures
                 * printed out.
                 * Think about ways around this in the future?
                 * For now, keep this choice in comments.
                 */
                /*bPres = (ir->eI==eiVV || IR_NPT_TROTTER(ir)); */
                /*bTemp = ((ir->eI==eiVV &&(!bInitStep)) || (ir->eI==eiVVAK && IR_NPT_TROTTER(ir)));*/
                bPres = TRUE;
                bTemp = ((ir->eI == eiVV && (!bInitStep)) || (ir->eI == eiVVAK));
                if (bCalcEner && ir->eI == eiVVAK)  /*MRS:  7/9/2010 -- this still doesn't fix it?*/
                {
                    bSumEkinhOld = TRUE;
                }
                /* for vv, the first half of the integration actually corresponds to the previous step.
                   So we need information from the last step in the first half of the integration */
                if (bGStat || do_per_step(step-1, nstglobalcomm))
                {
                    wallcycle_stop(wcycle, ewcUPDATE);
                    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                    wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                    constr, NULL, FALSE, state->box,
                                    top_global, &bSumEkinhOld,
                                    cglo_flags
                                    | CGLO_ENERGY
                                    | (bTemp ? CGLO_TEMPERATURE : 0)
                                    | (bPres ? CGLO_PRESSURE : 0)
                                    | (bPres ? CGLO_CONSTRAINT : 0)
                                    | (bStopCM ? CGLO_STOPCM : 0)
                                    | ((iterate.bIterationActive) ? CGLO_ITERATE : 0)
                                    | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                                    | CGLO_SCALEEKIN
                                    );
                    /* explanation of above:
                       a) We compute Ekin at the full time step
                       if 1) we are using the AveVel Ekin, and it's not the
                       initial step, or 2) if we are using AveEkin, but need the full
                       time step kinetic energy for the pressure (always true now, since we want accurate statistics).
                       b) If we are using EkinAveEkin for the kinetic energy for the temperature control, we still feed in
                       EkinAveVel because it's needed for the pressure */
                    wallcycle_start(wcycle, ewcUPDATE);
                }
                /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
                if (!bInitStep)
                {
                    if (bTrotter)
                    {
                        m_add(force_vir, shake_vir, total_vir); /* we need the un-dispersion corrected total vir here */
                        trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ2);
                    }
                }

                if (iterate.bIterationActive &&
                    done_iterating(cr, fplog, step, &iterate, bFirstIterate,
                                   state->veta, &vetanew))
                {
                    break;
                }
                bFirstIterate = FALSE;
            }

            if (bTrotter && !bInitStep)
            {
                copy_mat(shake_vir, state->svir_prev);
                copy_mat(force_vir, state->fvir_prev);
                if (IR_NVT_TROTTER(ir) && ir->eI == eiVV)
                {
                    /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
                    enerd->term[F_TEMP] = sum_ekin(&(ir->opts), ekind, NULL, (ir->eI == eiVV), FALSE);
                    enerd->term[F_EKIN] = trace(ekind->ekin);
                }
            }
            /* if it's the initial step, we performed this first step just to get the constraint virial */
            if (ir->eI == eiVV && bInitStep)
            {
                copy_rvecn(vbuf, state->v, 0, state->natoms);
                sfree(vbuf);
            }
            wallcycle_stop(wcycle, ewcUPDATE);
        }

        /* MRS -- now done iterating -- compute the conserved quantity */
        if (bVV)
        {
            saved_conserved_quantity = compute_conserved_from_auxiliary(ir, state, &MassQ);
            if (ir->eI == eiVV)
            {
                last_ekin = enerd->term[F_EKIN];
            }
            if ((ir->eDispCorr != edispcEnerPres) && (ir->eDispCorr != edispcAllEnerPres))
            {
                saved_conserved_quantity -= enerd->term[F_DISPCORR];
            }
        }

        /* ########  END FIRST UPDATE STEP  ############## */
        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         */
        do_md_trajectory_writing(fplog, cr, nfile, fnm, step, step_rel, t,
                                 ir, state, state_global, top_global, fr,
                                 outf, mdebin, ekind, f, f_global,
                                 &nchkpt,
                                 bCPT, FALSE, bLastStep, (Flags & MD_CONFOUT),
                                 bSumEkinhOld);

        /* Check if IMD step and do IMD communication, if bIMD is TRUE. */
        bIMDstep = do_IMD(ir->bIMD, step, cr, bNS, state->box, state->x, ir, t, wcycle);

        /* kludge -- virial is lost with restart for NPT control. Must restart */
        if (bStartingFromCpt && bVV)
        {
            copy_mat(state->svir_prev, shake_vir);
            copy_mat(state->fvir_prev, force_vir);
        }

        elapsed_time = walltime_accounting_get_current_elapsed_time(walltime_accounting);

        /* Check whether everything is still allright */
        if (((int)gmx_get_stop_condition() > handled_stop_condition)
#ifdef GMX_THREAD_MPI
            && MASTER(cr)
#endif
            )
        {
            /* this is just make gs.sig compatible with the hack
               of sending signals around by MPI_Reduce with together with
               other floats */
            if (gmx_get_stop_condition() == gmx_stop_cond_next_ns)
            {
                gs.sig[eglsSTOPCOND] = 1;
            }
            if (gmx_get_stop_condition() == gmx_stop_cond_next)
            {
                gs.sig[eglsSTOPCOND] = -1;
            }
            /* < 0 means stop at next step, > 0 means stop at next NS step */
            if (fplog)
            {
                fprintf(fplog,
                        "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                        gmx_get_signal_name(),
                        gs.sig[eglsSTOPCOND] == 1 ? "NS " : "");
                fflush(fplog);
            }
            fprintf(stderr,
                    "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                    gmx_get_signal_name(),
                    gs.sig[eglsSTOPCOND] == 1 ? "NS " : "");
            fflush(stderr);
            handled_stop_condition = (int)gmx_get_stop_condition();
        }
        else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
                 (max_hours > 0 && elapsed_time > max_hours*60.0*60.0*0.99) &&
                 gs.sig[eglsSTOPCOND] == 0 && gs.set[eglsSTOPCOND] == 0)
        {
            /* Signal to terminate the run */
            gs.sig[eglsSTOPCOND] = 1;
            if (fplog)
            {
                fprintf(fplog, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(step, sbuf), max_hours*0.99);
            }
            fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(step, sbuf), max_hours*0.99);
        }

        if (bResetCountersHalfMaxH && MASTER(cr) &&
            elapsed_time > max_hours*60.0*60.0*0.495)
        {
            gs.sig[eglsRESETCOUNTERS] = 1;
        }

        if (ir->nstlist == -1)
        {
            /* When bGStatEveryStep=FALSE, global_stat is only called
             * when we check the atom displacements, not at NS steps.
             * This means that also the bonded interaction count check is not
             * performed immediately after NS. Therefore a few MD steps could
             * be performed with missing interactions.
             * But wrong energies are never written to file,
             * since energies are only written after global_stat
             * has been called.
             */
            if (step >= nlh.step_nscheck)
            {
                nlh.nabnsb = natoms_beyond_ns_buffer(ir, fr, &top->cgs,
                                                     nlh.scale_tot, state->x);
            }
            else
            {
                /* This is not necessarily true,
                 * but step_nscheck is determined quite conservatively.
                 */
                nlh.nabnsb = 0;
            }
        }

        /* In parallel we only have to check for checkpointing in steps
         * where we do global communication,
         *  otherwise the other nodes don't know.
         */
        if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                           cpt_period >= 0 &&
                           (cpt_period == 0 ||
                            elapsed_time >= nchkpt*cpt_period*60.0)) &&
            gs.set[eglsCHKPT] == 0)
        {
            gs.sig[eglsCHKPT] = 1;
        }

        /* at the start of step, randomize or scale the velocities (trotter done elsewhere) */
        if (EI_VV(ir->eI))
        {
            if (!bInitStep)
            {
                update_tcouple(step - 1, ir, state, ekind, &MassQ, mdatoms);
            }
            if (ETC_ANDERSEN(ir->etc)) /* keep this outside of update_tcouple because of the extra info required to pass */
            {
                gmx_bool bIfRandomize;
                bIfRandomize = update_randomize_velocities(ir, step, cr, mdatoms, state, upd, constr);
                /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
                if (constr && bIfRandomize)
                {
                    update_constraints(fplog, step, NULL, ir, ekind, mdatoms,
                                       state, fr->bMolPBC, graph, f,
                                       &top->idef, tmp_vir,
                                       cr, nrnb, wcycle, upd, constr,
                                       TRUE, bCalcVir, vetanew);
                }
            }
        }

        if (bIterativeCase && do_per_step(step, ir->nstpcouple))
        {
            gmx_iterate_init(&iterate, TRUE);
            /* for iterations, we save these vectors, as we will be redoing the calculations */
            copy_coupling_state(state, bufstate, ekind, ekind_save, &(ir->opts));
        }

        /* for velocity Verlet integrator, we push at the end of a full step */
        if ( go->cfg->dohmc && EI_VV(ir->eI) ) {
          gmxgo_hmcpushv(go, state, step);
          gmxgo_hmcselect(go, state, f, 1, fr->ePBC, state->box, step);
        }

        bFirstIterate = TRUE;
        while (bFirstIterate || iterate.bIterationActive)
        {
            /* We now restore these vectors to redo the calculation with improved extended variables */
            if (iterate.bIterationActive)
            {
                copy_coupling_state(bufstate, state, ekind_save, ekind, &(ir->opts));
            }

            /* We make the decision to break or not -after- the calculation of Ekin and Pressure,
               so scroll down for that logic */

            /* #########   START SECOND UPDATE STEP ################# */
            /* Box is changed in update() when we do pressure coupling,
             * but we should still use the old box for energy corrections and when
             * writing it to the energy file, so it matches the trajectory files for
             * the same timestep above. Make a copy in a separate array.
             */
            copy_mat(state->box, lastbox);

            bOK         = TRUE;
            dvdl_constr = 0;

            {
                wallcycle_start(wcycle, ewcUPDATE);
                /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
                if (bTrotter)
                {
                    if (iterate.bIterationActive)
                    {
                        if (bFirstIterate)
                        {
                            scalevir = 1;
                        }
                        else
                        {
                            /* we use a new value of scalevir to converge the iterations faster */
                            scalevir = tracevir/trace(shake_vir);
                        }
                        msmul(shake_vir, scalevir, shake_vir);
                        m_add(force_vir, shake_vir, total_vir);
                        clear_mat(shake_vir);
                    }
                    trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ3);
                    /* We can only do Berendsen coupling after we have summed
                     * the kinetic energy or virial. Since the happens
                     * in global_state after update, we should only do it at
                     * step % nstlist = 1 with bGStatEveryStep=FALSE.
                     */
                }
                else
                {
                    update_tcouple(step, ir, state, ekind, &MassQ, mdatoms);
                    update_pcouple(fplog, step, ir, state, pcoupl_mu, M, bInitStep);
                }

                if (bVV)
                {
                    bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                    /* velocity half-step update */
                    update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                                  bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                                  ekind, M, upd, FALSE, etrtVELOCITY2,
                                  cr, nrnb, constr, &top->idef);
                }

                /* Above, initialize just copies ekinh into ekin,
                 * it doesn't copy position (for VV),
                 * and entire integrator for MD.
                 */

                if (ir->eI == eiVVAK)
                {
                    /* We probably only need md->homenr, not state->natoms */
                    if (state->natoms > cbuf_nalloc)
                    {
                        cbuf_nalloc = state->natoms;
                        srenew(cbuf, cbuf_nalloc);
                    }
                    copy_rvecn(state->x, cbuf, 0, state->natoms);
                }
                bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                              bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                              ekind, M, upd, bInitStep, etrtPOSITION, cr, nrnb, constr, &top->idef);
                wallcycle_stop(wcycle, ewcUPDATE);

                update_constraints(fplog, step, &dvdl_constr, ir, ekind, mdatoms, state,
                                   fr->bMolPBC, graph, f,
                                   &top->idef, shake_vir,
                                   cr, nrnb, wcycle, upd, constr,
                                   FALSE, bCalcVir, state->veta);

                /* we push the velocity here because of the following reason
                 * for the leapfrog integrator, the flow is shown below
                 *
                 *  1. neighbor search and domain decomposition
                 *  2. compute force          [ in do_force() ]
                 *    A. position x, and force f are pushed here for HMC
                 *
                 *  3. compute the velocity rescaling factor [ in update_tcouple() ]
                 *  4. v' = v + (f/m) dt      [ in update_coords() ]
                 *  5. x' = x + v dt          [ in update_coords() ]
                 *  6. constraints for x, v   [ in update_constraints ]
                 *    B. push v for HMC
                 *    C. possible HMC pop
                 *
                 * If an HMC pop occurs, it will follow through steps 1-2
                 * for the next MD step (which should not change the value of f).
                 * Next it should experience step 3, which requires the velocity
                 * to satisfy the constraints (eliminating normal components,
                 * which reduces the kinetic energy)
                 * Thus, a velocity push can only occur after step 2 or 6.
                 * In the latter case, as we proceed to step 4, we actually
                 * recover the pushed v' back to v, which fits natually into
                 * the velocity rescaling scheme, as shown below
                 *
                 *  --(step a2)--> x1,  v1    [push x1]
                 *  --(step a4)--> x1,  v2*
                 *  --(step a5)--> x2*, v2*
                 *  --(step a6)--> x2,  v2    [push v2]
                 * ---------------------------------
                 *   normal MD propagation
                 * ---------------------------------
                 *  --(step b2)--> x2,  v2    [push x2]
                 *  --(step b4)--> x2,  v3*
                 *  --(step b5)--> x3*, v3*
                 *  --(step b6)--> x3,  v3    [push v3]
                 * ----------------------------------
                 *   HMC pop x2, -v3, f2
                 *   roughly recovers the state after b4
                 *   with the inverted and normalized v3
                 * ----------------------------------
                 *  --(step c2)--> x2, -v3
                 *  --(step c4)--> x2, -v2#   [reversing b4, roughly the state after b2]
                 *  --(step c5)--> x1#,-v2#   [reversing a5, roughly the state after a4, except the force]
                 *  --(step c6)--> x1, -v2
                 *
                 * For the forward propagation
                 *  In step a5,
                 *    x2* = x1 + (v2*) dt     [x.forward] <-------------+
                 *  In step b4,                                         |
                 *    v3* = v2 + (f2/m) dt    [v.forward] <--------+    |
                 *  In step b6,                                    |    |
                 *    v3 = Normalize(v3*)                          |    |
                 *                                                 |    |
                 * For the reverse,                                |    |
                 *  In step c4,                                    |    |
                 *   -v2# = -v3 + (f2/m) dt   [v.backward] <-------+    |
                 *  In step c5,                                         |
                 *    x1# = x2 + (-v2#) dt    [x.backward] <------------+
                 *  In step c6,
                 *   -v2  = Normalize(-v2#)
                 *
                 * Another strategy is to place the push before constraints
                 * and manually disable the velocity rescaling (thermostat)
                 * in the next step.  But this seems to be an overkill,
                 * unless the current scheme fails.
                 * */
                if ( go->cfg->dohmc && !EI_VV(ir->eI) ) {
                  gmxgo_hmcpushv(go, state, step);
                }

                if (bCalcVir && bUpdateDoLR && ir->nstcalclr > 1)
                {
                    /* Correct the virial for multiple time stepping */
                    m_sub(shake_vir, fr->vir_twin_constr, shake_vir);
                }

                if (ir->eI == eiVVAK)
                {
                    /* erase F_EKIN and F_TEMP here? */
                    /* just compute the kinetic energy at the half step to perform a trotter step */
                    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                    wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                    constr, NULL, FALSE, lastbox,
                                    top_global, &bSumEkinhOld,
                                    cglo_flags | CGLO_TEMPERATURE
                                    );
                    wallcycle_start(wcycle, ewcUPDATE);
                    trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ4);
                    /* now we know the scaling, we can compute the positions again again */
                    copy_rvecn(cbuf, state->x, 0, state->natoms);

                    bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                    update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                                  bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                                  ekind, M, upd, bInitStep, etrtPOSITION, cr, nrnb, constr, &top->idef);
                    wallcycle_stop(wcycle, ewcUPDATE);

                    /* do we need an extra constraint here? just need to copy out of state->v to upd->xp? */
                    /* are the small terms in the shake_vir here due
                     * to numerical errors, or are they important
                     * physically? I'm thinking they are just errors, but not completely sure.
                     * For now, will call without actually constraining, constr=NULL*/
                    update_constraints(fplog, step, NULL, ir, ekind, mdatoms,
                                       state, fr->bMolPBC, graph, f,
                                       &top->idef, tmp_vir,
                                       cr, nrnb, wcycle, upd, NULL,
                                       FALSE, bCalcVir,
                                       state->veta);
                }
                if (!bOK)
                {
                    gmx_fatal(FARGS, "Constraint error: Shake, Lincs or Settle could not solve the constraints");
                }

                if (fr->bSepDVDL && fplog && do_log)
                {
                    gmx_print_sepdvdl(fplog, "Constraint dV/dl", 0.0, dvdl_constr);
                }
                if (bVV)
                {
                    /* this factor or 2 correction is necessary
                       because half of the constraint force is removed
                       in the vv step, so we have to double it.  See
                       the Redmine issue #1255.  It is not yet clear
                       if the factor of 2 is exact, or just a very
                       good approximation, and this will be
                       investigated.  The next step is to see if this
                       can be done adding a dhdl contribution from the
                       rattle step, but this is somewhat more
                       complicated with the current code. Will be
                       investigated, hopefully for 4.6.3. However,
                       this current solution is much better than
                       having it completely wrong.
                     */
                    enerd->term[F_DVDL_CONSTR] += 2*dvdl_constr;
                }
                else
                {
                    enerd->term[F_DVDL_CONSTR] += dvdl_constr;
                }
            }

            /* ############## IF NOT VV, Calculate globals HERE, also iterate constraints  ############ */
            /* With Leap-Frog we can skip compute_globals at
             * non-communication steps, but we need to calculate
             * the kinetic energy one step before communication.
             */
            if (bGStat || (!EI_VV(ir->eI) && do_per_step(step+1, nstglobalcomm)))
            {
                if (ir->nstlist == -1 && bFirstIterate)
                {
                    gs.sig[eglsNABNSB] = nlh.nabnsb;
                }
                compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                constr,
                                bFirstIterate ? &gs : NULL,
                                (step_rel % gs.nstms == 0) &&
                                (multisim_nsteps < 0 || (step_rel < multisim_nsteps)),
                                lastbox,
                                top_global, &bSumEkinhOld,
                                cglo_flags
                                | (!EI_VV(ir->eI) ? CGLO_ENERGY : 0)
                                | (!EI_VV(ir->eI) && bStopCM ? CGLO_STOPCM : 0)
                                | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                                | (!EI_VV(ir->eI) ? CGLO_PRESSURE : 0)
                                | (iterate.bIterationActive ? CGLO_ITERATE : 0)
                                | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                                | CGLO_CONSTRAINT
                                );
                if (ir->nstlist == -1 && bFirstIterate)
                {
                    nlh.nabnsb         = gs.set[eglsNABNSB];
                    gs.set[eglsNABNSB] = 0;
                }
            }
            /* bIterate is set to keep it from eliminating the old ekin kinetic energy terms */
            /* #############  END CALC EKIN AND PRESSURE ################# */

            /* Note: this is OK, but there are some numerical precision issues with using the convergence of
               the virial that should probably be addressed eventually. state->veta has better properies,
               but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
               generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

            if (iterate.bIterationActive &&
                done_iterating(cr, fplog, step, &iterate, bFirstIterate,
                               trace(shake_vir), &tracevir))
            {
                break;
            }
            bFirstIterate = FALSE;
        }

        if (!bVV)
        {
            /* sum up the foreign energy and dhdl terms for md and sd. currently done every step so that dhdl is correct in the .edr */
            sum_dhdl(enerd, state->lambda, ir->fepvals);
        }
        update_box(fplog, step, ir, mdatoms, state, f,
                   ir->nstlist == -1 ? &nlh.scale_tot : NULL, pcoupl_mu, nrnb, upd);

        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */

        if ( go->cfg->dohmc && !EI_VV(ir->eI) ) {
          gmxgo_hmcselect(go, state, f, 1, fr->ePBC, state->box, step);
        }

        /* The coordinates (x) were unshifted in update */
        if (!bGStat)
        {
            /* We will not sum ekinh_old,
             * so signal that we still have to do it.
             */
            bSumEkinhOld = TRUE;
        }

        /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

        /* use the directly determined last velocity, not actually the averaged half steps */
        if (bTrotter && ir->eI == eiVV)
        {
            enerd->term[F_EKIN] = last_ekin;
        }
        enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

        if (bVV)
        {
            enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + saved_conserved_quantity;
        }
        else
        {
            enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + compute_conserved_from_auxiliary(ir, state, &MassQ);
        }
        /* #########  END PREPARING EDR OUTPUT  ###########  */

        /* Output stuff */
        if (MASTER(cr))
        {
            gmx_bool do_dr, do_or;

            if (!(bStartingFromCpt && (EI_VV(ir->eI))))
            {
                if (bCalcEner)
                {
                    upd_mdebin(mdebin, FALSE, TRUE,
                               t, mdatoms->tmass, enerd, state,
                               ir->fepvals, ir->expandedvals, lastbox,
                               shake_vir, force_vir, total_vir, pres,
                               ekind, mu_tot, constr);
                }
                else
                {
                    upd_mdebin_step(mdebin);
                }

                do_dr  = do_per_step(step, ir->nstdisreout);
                do_or  = do_per_step(step, ir->nstorireout);

                print_ebin(mdoutf_get_fp_ene(outf), do_ene, do_dr, do_or, do_log ? fplog : NULL,
                           step, t,
                           eprNORMAL, bCompact, mdebin, fcd, groups, &(ir->opts));
            }

            if (do_per_step(step, ir->nstlog))
            {
                if (fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
                }
            }
        }
        /* Print the remaining wall clock time for the run */
        if (MULTIMASTER(cr) && (do_verbose || gmx_got_usr_signal()) && !bPMETuneRunning)
        {
            print_time(stderr, walltime_accounting, step, ir, cr);
        }

        /* Ion/water position swapping.
         * Not done in last step since trajectory writing happens before this call
         * in the MD loop and exchanges would be lost anyway. */
        bNeedRepartition = FALSE;
        if ((ir->eSwapCoords != eswapNO) && (step > 0) && !bLastStep &&
            do_per_step(step, ir->swap->nstswap))
        {
            bNeedRepartition = do_swapcoords(cr, step, t, ir, wcycle,
                                             state->x,
                                             state->box,
                                             top_global, MASTER(cr) && bVerbose, FALSE);

            if (bNeedRepartition && DOMAINDECOMP(cr))
            {
                dd_collect_state(cr->dd, state, state_global);
            }
        }

        if ( bNeedRepartition && DOMAINDECOMP(cr) )
        {
            dd_partition_system(fplog, step, cr, TRUE, 1,
                                state_global, top_global, ir,
                                state, &f, mdatoms, top, fr,
                                NULL, NULL, constr,
                                nrnb, wcycle, FALSE);
        }

        bFirstStep       = FALSE;
        bInitStep        = FALSE;
        bStartingFromCpt = FALSE;

        /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
        /* With all integrators, except VV, we need to retain the pressure
         * at the current step for coupling at the next step.
         */
        if ((state->flags & (1<<estPRES_PREV)) &&
            (bGStatEveryStep ||
             (ir->nstpcouple > 0 && step % ir->nstpcouple == 0)))
        {
            /* Store the pressure in t_state for pressure coupling
             * at the next MD step.
             */
            copy_mat(pres, state->pres_prev);
        }

        /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */

        {
            /* increase the MD step number */
            step++;
            step_rel++;
        }

        cycles = wallcycle_stop(wcycle, ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        if (bPMETuneRunning || bPMETuneTry)
        {
            /* PME grid + cut-off optimization with GPUs or PME nodes */

            /* Count the total cycles over the last steps */
            cycles_pmes += cycles;

            /* We can only switch cut-off at NS steps */
            if (step % ir->nstlist == 0)
            {
                /* PME grid + cut-off optimization with GPUs or PME nodes */
                if (bPMETuneTry)
                {
                    if (DDMASTER(cr->dd))
                    {
                        /* PME node load is too high, start tuning */
                        bPMETuneRunning = (dd_pme_f_ratio(cr->dd) >= 1.05);
                    }
                    dd_bcast(cr->dd, sizeof(gmx_bool), &bPMETuneRunning);

                    if (bPMETuneRunning &&
                        fr->nbv->bUseGPU && DOMAINDECOMP(cr) &&
                        !(cr->duty & DUTY_PME))
                    {
                        /* Lock DLB=auto to off (does nothing when DLB=yes/no).
                         * With GPUs + separate PME ranks, we don't want DLB.
                         * This could happen when we scan coarse grids and
                         * it would then never be turned off again.
                         * This would hurt performance at the final, optimal
                         * grid spacing, where DLB almost never helps.
                         * Also, DLB can limit the cut-off for PME tuning.
                         */
                        dd_dlb_set_lock(cr->dd, TRUE);
                    }

                    if (bPMETuneRunning || step_rel > ir->nstlist*50)
                    {
                        bPMETuneTry     = FALSE;
                    }
                }
                if (bPMETuneRunning)
                {
                    /* init_step might not be a multiple of nstlist,
                     * but the first cycle is always skipped anyhow.
                     */
                    bPMETuneRunning =
                        pme_load_balance(pme_loadbal, cr,
                                         (bVerbose && MASTER(cr)) ? stderr : NULL,
                                         fplog,
                                         ir, state, cycles_pmes,
                                         fr->ic, fr->nbv, &fr->pmedata,
                                         step);

                    /* Update constants in forcerec/inputrec to keep them in sync with fr->ic */
                    fr->ewaldcoeff_q  = fr->ic->ewaldcoeff_q;
                    fr->ewaldcoeff_lj = fr->ic->ewaldcoeff_lj;
                    fr->rlist         = fr->ic->rlist;
                    fr->rlistlong     = fr->ic->rlistlong;
                    fr->rcoulomb      = fr->ic->rcoulomb;
                    fr->rvdw          = fr->ic->rvdw;

                    if (ir->eDispCorr != edispcNO)
                    {
                        calc_enervirdiff(NULL, ir->eDispCorr, fr);
                    }

                    if (!bPMETuneRunning &&
                        DOMAINDECOMP(cr) &&
                        dd_dlb_is_locked(cr->dd))
                    {
                        /* Unlock the DLB=auto, DLB is allowed to activate
                         * (but we don't expect it to activate in most cases).
                         */
                        dd_dlb_set_lock(cr->dd, FALSE);
                    }
                }
                cycles_pmes = 0;
            }
        }

        if (step_rel == wcycle_get_reset_counters(wcycle) ||
            gs.set[eglsRESETCOUNTERS] != 0)
        {
            /* Reset all the counters related to performance over the run */
            reset_all_counters(fplog, cr, step, &step_rel, ir, wcycle, nrnb, walltime_accounting,
                               fr->nbv != NULL && fr->nbv->bUseGPU ? fr->nbv->cu_nbv : NULL);
            wcycle_set_reset_counters(wcycle, -1);
            if (!(cr->duty & DUTY_PME))
            {
                /* Tell our PME node to reset its counters */
                gmx_pme_send_resetcounters(cr, step);
            }
            /* Correct max_hours for the elapsed time */
            max_hours                -= elapsed_time/(60.0*60.0);
            bResetCountersHalfMaxH    = FALSE;
            gs.set[eglsRESETCOUNTERS] = 0;
        }

        /* If bIMD is TRUE, the master updates the IMD energy record and sends positions to VMD client */
        IMD_prep_energies_send_positions(ir->bIMD && MASTER(cr), bIMDstep, ir->imd, enerd, step, bCalcEner, wcycle);

    }
    /* End of main MD loop */
    debug_gmx();

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end(walltime_accounting);

    if (!(cr->duty & DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0)
        {
            print_ebin(mdoutf_get_fp_ene(outf), FALSE, FALSE, FALSE, fplog, step, t,
                       eprAVER, FALSE, mdebin, fcd, groups, &(ir->opts));
        }
    }

    done_mdoutf(outf);
    debug_gmx();

    if (ir->nstlist == -1 && nlh.nns > 0 && fplog)
    {
        fprintf(fplog, "Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n", nlh.s1/nlh.nns, sqrt(nlh.s2/nlh.nns - sqr(nlh.s1/nlh.nns)));
        fprintf(fplog, "Average number of atoms that crossed the half buffer length: %.1f\n\n", nlh.ab/nlh.nns);
    }

    if (pme_loadbal != NULL)
    {
        pme_loadbal_done(pme_loadbal, cr, fplog,
                         fr->nbv != NULL && fr->nbv->bUseGPU);
    }

    /* IMD cleanup, if bIMD is TRUE. */
    IMD_finalize(ir->bIMD, ir->imd);

    walltime_accounting_set_nsteps_done(walltime_accounting, step_rel);

    gmxgo_close(go);

    return 0;
}
