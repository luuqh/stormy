! $Id: climat.h 697 2011-04-11 12:35:17Z gcambon $
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://roms.mpl.ird.fr
!======================================================================
!
/*  This is include file "climat.h"
-----------------------------------
 Free surface climatology:
 ==== ======= ============
   ssh        sea surface height climatology at current time-step.
   Znudgcof   inverse relaxation time [1/sec] for nudging toward
                               free surface climatological fields.
   sshg       two-time-level array to hold climatological data for 
                                                     free surface.
   tssh       time of read in sea surface height climatology.
*/
#if defined ZCLIMATOLOGY || defined AGRIF
      real ssh(GLOBAL_2D_ARRAY)
      common /climat_ssh/ssh
#endif
#ifdef ZCLIMATOLOGY
# ifdef ZNUDGING
      real Znudgcof(GLOBAL_2D_ARRAY)
      common /climat_Znudgcof/Znudgcof
# endif
# ifndef ANA_SSH
      real sshg(GLOBAL_2D_ARRAY,2)
      common /climat_sshg/sshg

      real    ssh_time(2)
      real    ssh_cycle
      integer itssh, ssh_ncycle, ssh_rec, ssh_tid, ssh_id 
      common /climat_zdat1/ ssh_time
      common /climat_zdat2/ ssh_cycle
      common /climat_zdat3/ 
     &        itssh, ssh_ncycle, ssh_rec, ssh_tid, ssh_id

#   undef SSH_DATA
# endif /* !ANA_SSH */
#endif

/*
 Temperature and salinity climatology:
 =========== === ======== ============
   tclm       climatology for tracer variables at current time-step.
   Tnudgcof   inverse relaxation time [1/sec] for nudging toward
                                       tracer climatological fields.
   tclima     two-time-level array to hold climatological data for
                                               tracer variables.
   ttclm      time of read in climatology for tracer type variables.
*/
#ifdef SOLVE3D
# if defined TCLIMATOLOGY || (defined AGRIF && !defined T_FRC_BRY)
      real tclm(GLOBAL_2D_ARRAY,N,NT)
      common /climat_tclm/tclm
# endif
# ifdef TCLIMATOLOGY
#  ifdef TNUDGING
      real Tnudgcof(GLOBAL_2D_ARRAY,N,NT)
      common /climat_Tnudgcof/Tnudgcof
#  endif
#  ifndef ANA_TCLIMA
      real tclima(GLOBAL_2D_ARRAY,N,2,NT)
      common /climat_tclima/tclima

      real tclm_time(2,NT)
      real tclm_cycle(NT)
      integer ittclm(NT), tclm_ncycle(NT), tclm_rec(NT), 
     &        tclm_tid(NT), tclm_id(NT)
      logical got_tclm(NT)

      common /climat_tdat/  tclm_time,       tclm_cycle,
     &        ittclm,       tclm_ncycle,     tclm_rec,
     &                      tclm_tid,        tclm_id,
     &                                       got_tclm

#   undef TCLIMA_DATA
#  endif /* !ANA_TCLIMA */
# endif /* TCLIMATOLOGY */
#endif /* SOLVE3D */
/*
 barotropic and baroclinic velocity climatology:
 ========== === ========== ======== ===========
   ubclm     climatology for bar. u-velocity at current time-step.
   vbclm     climatology for bar. v-velocity at current time-step.
   uclm      climatology for u-velocity at current time-step.
   vclm      climatology for v-velocity at current time-step.

   ubclima   two-time-level array to hold climatological data
   vbclima
   uclima
   vclima
*/
#if defined M2CLIMATOLOGY || (defined AGRIF && !defined M2_FRC_BRY)
      real ubclm(GLOBAL_2D_ARRAY)
      real vbclm(GLOBAL_2D_ARRAY)
      common /climat_ubclm/ubclm /climat_vbclm/vbclm 
#endif
#if defined SOLVE3D && (defined M3CLIMATOLOGY || \
                        (defined AGRIF && !defined M3_FRC_BRY))
      real uclm(GLOBAL_2D_ARRAY,N)
      real vclm(GLOBAL_2D_ARRAY,N)
      common /climat_uclm/uclm /climat_vclm/vclm
#endif
#ifdef M2CLIMATOLOGY
# ifdef M2NUDGING
      real M2nudgcof(GLOBAL_2D_ARRAY)
      common /climat_M2nudgcof/M2nudgcof
# endif
# ifndef ANA_M2CLIMA
      real ubclima(GLOBAL_2D_ARRAY,2)
      real vbclima(GLOBAL_2D_ARRAY,2)
      common /climat_ubclima/ubclima /climat_vbclima/vbclima
# endif
#endif
!
#if defined SOLVE3D && defined M3CLIMATOLOGY
#   ifdef M3NUDGING
      real M3nudgcof(GLOBAL_2D_ARRAY)
      common /climat_M3nudgcof/M3nudgcof
#   endif
#   ifndef ANA_M3CLIMA
      real uclima(GLOBAL_2D_ARRAY,N,2)
      real vclima(GLOBAL_2D_ARRAY,N,2)
      common /climat_uclima/uclima /climat_vclima/vclima
#   endif
#endif
!
#if defined M2CLIMATOLOGY || defined M3CLIMATOLOGY
      real     uclm_time(2)
      real     uclm_cycle
      integer ituclm, uclm_ncycle, uclm_rec, uclm_tid,
     &        ubclm_id, vbclm_id, uclm_id, vclm_id
      common /climat_udat1/  uclm_time
      common /climat_udat2/  uclm_cycle
      common /climat_udat3/
     &             ituclm,   uclm_ncycle, uclm_rec,
     &             uclm_tid, ubclm_id,    vbclm_id,
     &             uclm_id,  vclm_id
!
#endif

