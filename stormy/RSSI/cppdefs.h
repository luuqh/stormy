! $Id: cppdefs.h 773 2012-02-01 17:05:30Z marchesiello $
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
/*
   This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ====
*/
#undef  BASIN           /* Basin Example */
#undef  CANYON_A        /* Canyon_A Example */
#undef  CANYON_B        /* Canyon_B Example */
#undef  EQUATOR         /* Equator Example  */
#undef  GRAV_ADJ        /* Graviational Adjustment Example */
#undef  INNERSHELF      /* Inner Shelf Example */
#undef  RIVER           /* River run-off Example */
#undef  OVERFLOW        /* Graviational/Overflow Example */
#undef  SEAMOUNT        /* Seamount Example */
#undef  SHELFRONT       /* Shelf Front Example */
#undef  SOLITON         /* Equatorial Rossby Wave Example */
#undef  UPWELLING       /* Upwelling Example */
#undef  VORTEX          /* Baroclinic Vortex Example */
#undef  INTERNAL        /* Internal Tide Example */
#define REGIONAL        /* REGIONAL Applications */


#if defined REGIONAL
/*
!====================================================================
!               REGIONAL (realistic) Configurations
!==================================================================== 
!
!------------------------
! BASIC OPTIONS
!------------------------
!
*/
                      /* Configuration Name */
# define STSU
                      /* Parallelization */
# undef  OPENMP
# define  MPI
                      /* Nesting */
# undef  AGRIF
# undef  AGRIF_2WAY
                      /* Open Boundary Conditions */
# undef  TIDES
# define OBC_EAST
# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH
                      /* Applications */
# undef  BIOLOGY
# undef  FLOATS
# undef  STATIONS
# undef  PASSIVE_TRACER
# undef  SEDIMENT
# undef  BBL
/*!
!------------------------
! PRE-SELECTED OPTIONS
!------------------------
*/
                      /* Parallelization */
# ifdef MPI
#  undef  PARALLEL_FILES
# endif
# undef  AUTOTILING
# undef  ETALON_CHECK
                      /* Model dynamics */
# undef SOLVE3D
# define UV_COR
# define UV_ADV
# ifdef TIDES
#  define SSH_TIDES
#  define UV_TIDES
#  define TIDERAMP
# endif
                      /* Grid configuration */
# define CURVGRID
# define SPHERICAL
# define MASKING
                      /* Lateral Momentum Mixing */
# define UV_VIS2
# define MIX_GP_UV
# undef  VIS_SMAGO
                      /* Lateral Tracer Mixing */
# undef MIX_GP_TS
# undef TS_DIF2
# undef  TS_SPLIT_UP3
                      /* Vertical Mixing */
# undef  BODYFORCE
# undef  BVF_MIXING
# undef LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_SKPP
#  define LMD_SKPP2005
#  define LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
#  undef  LMD_DDMIX
#  undef LMD_NONLOCAL
# endif
                      /* Equation of State */
# undef SALINITY
# undef NONLIN_EOS
# undef SPLIT_EOS
                      /* Surface Forcing */
# undef  BULK_FLUX
# ifdef BULK_FLUX
#  define BULK_FAIRALL
#  define BULK_LW
#  define BULK_EP
#  define BULK_SMFLUX
# else
#  undef QCORRECTION
#  undef SFLX_CORR
#  undef DIURNAL_SRFLUX
# endif
                      /* Lateral Forcing */
# define SPONGE

# undef CLIMATOLOGY
# ifdef CLIMATOLOGY
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define TCLIMATOLOGY

#  define ZNUDGING
#  define M2NUDGING
#  define M3NUDGING
#  define TNUDGING
#  undef  ROBUST_DIAG
# endif

# define  FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  undef M3_FRC_BRY
#  undef T_FRC_BRY
# endif
                      /* Bottom Forcing */
# define ANA_BSFLUX
# define ANA_BTFLUX
                      /* Storm Surge 2D model */
# undef ANA_BMFLUX
# define ANA_INITIAL
# define ANA_SRFLUX
# define ANA_BRY

# define ATM_PRESS

                      /* Point Sources - Rivers */
# undef  PSOURCE
# undef  ANA_PSOURCE
                      /* Open Boundary Conditions */
# ifdef TIDES
#  define OBC_M2FLATHER
# else
#  undef  OBC_M2SPECIFIED
#  undef  OBC_M2FLATHER
#  define OBC_M2CHARACT
#  undef  OBC_M2ORLANSKI
#  ifdef  OBC_M2ORLANSKI
#   define OBC_VOLCONS
#  endif
# endif
# undef OBC_M3ORLANSKI
# undef OBC_TORLANSKI
# undef  OBC_M3SPECIFIED
# undef  OBC_TSPECIFIED
                      /* Input/Output & Diagnostics */
# define AVERAGES
# undef AVERAGES_K
# undef  DIAGNOSTICS_TS
# undef  DIAGNOSTICS_UV
# ifdef DIAGNOSTICS_TS
#  undef DIAGNOSTICS_TS_ADV
#  undef DIAGNOSTICS_TS_MLD
# endif
/*
!           Applications:
!---------------------------------
! Biology, floats, Stations, 
! Passive tracer, Sediments, BBL
!---------------------------------
*/
                      /*      Choice of Biology models   */
# ifdef BIOLOGY
#  undef  PISCES
#  define BIO_NChlPZD
#  undef  BIO_N2P2Z2D2
#  undef  BIO_N2ChlPZD2  
                      /*  Options  */
#  ifdef PISCES
#   define key_trc_pisces
#   define key_passivetrc
#   undef  DIAGNOSTICS_BIO
#   ifdef DIAGNOSTICS_BIO
#     define key_trc_diaadd
#     define key_trc_dia3d
#   endif
#  endif
#  ifdef BIO_NChlPZD
#   undef  OXYGEN
#  endif
#  ifdef BIO_NChlPZD
#   define DIAGNOSTICS_BIO
#  endif
#  ifdef BIO_N2P2Z2D2
#   undef  VAR_CHL_C
#  endif
# endif
                      /*     Lagrangian floats model    */
# ifdef FLOATS
#  undef  FLOATS_GLOBAL_ATTRIBUTES
#  undef  IBM
#  undef  RANDOM_WALK
#  ifdef RANDOM_WALK
#   define DIEL_MIGRATION
#   define RANDOM_VERTICAL
#   define RANDOM_HORIZONTAL
#  endif
# endif
                      /*     Stations recording    */
# ifdef STATIONS
#  define ALL_SIGMA
# endif
                      /*      Sediment dynamics model     */
# ifdef SEDIMENT
#  define ANA_SEDIMENT
#  undef  BED_ARMOR
#  undef  ANA_SPFLUX
#  undef  ANA_BPFLUX
#  define LINEAR_CONTINUATION
#  undef  NEUMANN
# endif
                      /*      Bottom Boundary Layer model     */
# ifdef BBL
#  define ANA_WWAVE
#  ifdef SEDIMENT
#   undef  ANA_BSEDIM
#  else
#   define ANA_BSEDIM
#  endif
#  undef  Z0_BL
#  ifdef Z0_BL
#   define Z0_RIP
#  endif
#  undef  Z0_BIO
# endif
# endif
#include "set_global_definitions.h"

