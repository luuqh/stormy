#include "cppdefs.h"


      subroutine check_switches1 (ierr)

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if new CPP-switches  !!!!!!
!!!!!! were introduced into "cppdefs.h".  NO ACTION IS NEEDED  !!!!!!
!!!!!! if changes in "cppdefs.h" are limited to activation or  !!!!!!
!!!!!! deactivation of previously known switches.              !!!!!!
!!!!!! To refresh this file compile and execute "cppcheck.F"   !!!!!!
!!!!!! as an independent program, or use commands              !!!!!!
!!!!!! "make checkdefs" or "make depend".                      !!!!!!
!!!!!! Number of Configuration Choices:                     15 !!!!!!
!!!!!! Total number of CPP-switches:                       239 !!!!!!

      implicit none
      integer ierr, is,ie, iexample
#include "param.h"
#include "strings.h"
#ifdef MPI
# include "scalars.h"
#endif
      MPI_master_only write(stdout,'(/1x,A/)')
     &      'Activated C-preprocessing Options:'
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      iexample=0
      is=1
#ifdef BASIN
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'BASIN'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BASIN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CANYON_A
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'CANYON_A'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CANYON_A'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CANYON_B
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'CANYON_B'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CANYON_B'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef EQUATOR
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'EQUATOR'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='EQUATOR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GRAV_ADJ
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'GRAV_ADJ'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GRAV_ADJ'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef INNERSHELF
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'INNERSHELF'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='INNERSHELF'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RIVER
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'RIVER'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RIVER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OVERFLOW
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'OVERFLOW'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OVERFLOW'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SEAMOUNT
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'SEAMOUNT'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SEAMOUNT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SHELFRONT
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'SHELFRONT'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SHELFRONT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SOLITON
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'SOLITON'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SOLITON'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UPWELLING
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'UPWELLING'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UPWELLING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VORTEX
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'VORTEX'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VORTEX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef INTERNAL
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'INTERNAL'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='INTERNAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef REGIONAL
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'REGIONAL'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='REGIONAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef STSU
      MPI_master_only write(stdout,'(10x,A)') 'STSU'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='STSU'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OPENMP
      MPI_master_only write(stdout,'(10x,A)') 'OPENMP'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OPENMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MPI
      MPI_master_only write(stdout,'(10x,A)') 'MPI'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_2WAY
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_2WAY'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_2WAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TIDES
      MPI_master_only write(stdout,'(10x,A)') 'TIDES'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TIDES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_EAST
      MPI_master_only write(stdout,'(10x,A)') 'OBC_EAST'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_EAST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_WEST
      MPI_master_only write(stdout,'(10x,A)') 'OBC_WEST'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_WEST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_NORTH
      MPI_master_only write(stdout,'(10x,A)') 'OBC_NORTH'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_NORTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_SOUTH
      MPI_master_only write(stdout,'(10x,A)') 'OBC_SOUTH'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_SOUTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BIOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'BIOLOGY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BIOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef FLOATS
      MPI_master_only write(stdout,'(10x,A)') 'FLOATS'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FLOATS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef STATIONS
      MPI_master_only write(stdout,'(10x,A)') 'STATIONS'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='STATIONS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PASSIVE_TRACER
      MPI_master_only write(stdout,'(10x,A)') 'PASSIVE_TRACER'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PASSIVE_TRACER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SEDIMENT
      MPI_master_only write(stdout,'(10x,A)') 'SEDIMENT'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SEDIMENT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BBL
      MPI_master_only write(stdout,'(10x,A)') 'BBL'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BBL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PARALLEL_FILES
      MPI_master_only write(stdout,'(10x,A)') 'PARALLEL_FILES'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PARALLEL_FILES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AUTOTILING
      MPI_master_only write(stdout,'(10x,A)') 'AUTOTILING'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AUTOTILING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ETALON_CHECK
      MPI_master_only write(stdout,'(10x,A)') 'ETALON_CHECK'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ETALON_CHECK'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SOLVE3D
      MPI_master_only write(stdout,'(10x,A)') 'SOLVE3D'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SOLVE3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_COR
      MPI_master_only write(stdout,'(10x,A)') 'UV_COR'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_COR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_ADV
      MPI_master_only write(stdout,'(10x,A)') 'UV_ADV'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_ADV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SSH_TIDES
      MPI_master_only write(stdout,'(10x,A)') 'SSH_TIDES'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SSH_TIDES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_TIDES
      MPI_master_only write(stdout,'(10x,A)') 'UV_TIDES'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_TIDES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TIDERAMP
      MPI_master_only write(stdout,'(10x,A)') 'TIDERAMP'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TIDERAMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CURVGRID
      MPI_master_only write(stdout,'(10x,A)') 'CURVGRID'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CURVGRID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPHERICAL
      MPI_master_only write(stdout,'(10x,A)') 'SPHERICAL'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPHERICAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MASKING
      MPI_master_only write(stdout,'(10x,A)') 'MASKING'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MASKING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_VIS2
      MPI_master_only write(stdout,'(10x,A)') 'UV_VIS2'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VIS2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MIX_GP_UV
      MPI_master_only write(stdout,'(10x,A)') 'MIX_GP_UV'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MIX_GP_UV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VIS_SMAGO
      MPI_master_only write(stdout,'(10x,A)') 'VIS_SMAGO'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VIS_SMAGO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MIX_GP_TS
      MPI_master_only write(stdout,'(10x,A)') 'MIX_GP_TS'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MIX_GP_TS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_DIF2
      MPI_master_only write(stdout,'(10x,A)') 'TS_DIF2'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_DIF2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_SPLIT_UP3
      MPI_master_only write(stdout,'(10x,A)') 'TS_SPLIT_UP3'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_SPLIT_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BODYFORCE
      MPI_master_only write(stdout,'(10x,A)') 'BODYFORCE'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BODYFORCE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BVF_MIXING
      MPI_master_only write(stdout,'(10x,A)') 'BVF_MIXING'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BVF_MIXING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_MIXING
      MPI_master_only write(stdout,'(10x,A)') 'LMD_MIXING'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_MIXING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_SKPP
      MPI_master_only write(stdout,'(10x,A)') 'LMD_SKPP'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_SKPP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_SKPP2005
      MPI_master_only write(stdout,'(10x,A)') 'LMD_SKPP2005'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_SKPP2005'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_BKPP
      MPI_master_only write(stdout,'(10x,A)') 'LMD_BKPP'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_BKPP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_RIMIX
      MPI_master_only write(stdout,'(10x,A)') 'LMD_RIMIX'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_RIMIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_CONVEC
      MPI_master_only write(stdout,'(10x,A)') 'LMD_CONVEC'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_CONVEC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_DDMIX
      MPI_master_only write(stdout,'(10x,A)') 'LMD_DDMIX'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_DDMIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_NONLOCAL
      MPI_master_only write(stdout,'(10x,A)') 'LMD_NONLOCAL'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_NONLOCAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SALINITY
      MPI_master_only write(stdout,'(10x,A)') 'SALINITY'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SALINITY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NONLIN_EOS
      MPI_master_only write(stdout,'(10x,A)') 'NONLIN_EOS'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NONLIN_EOS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPLIT_EOS
      MPI_master_only write(stdout,'(10x,A)') 'SPLIT_EOS'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPLIT_EOS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_FLUX
      MPI_master_only write(stdout,'(10x,A)') 'BULK_FLUX'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_FLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_FAIRALL
      MPI_master_only write(stdout,'(10x,A)') 'BULK_FAIRALL'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_FAIRALL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_LW
      MPI_master_only write(stdout,'(10x,A)') 'BULK_LW'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_LW'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_EP
      MPI_master_only write(stdout,'(10x,A)') 'BULK_EP'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_EP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_SMFLUX
      MPI_master_only write(stdout,'(10x,A)') 'BULK_SMFLUX'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_SMFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef QCORRECTION
      MPI_master_only write(stdout,'(10x,A)') 'QCORRECTION'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='QCORRECTION'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SFLX_CORR
      MPI_master_only write(stdout,'(10x,A)') 'SFLX_CORR'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SFLX_CORR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIURNAL_SRFLUX
      MPI_master_only write(stdout,'(10x,A)') 'DIURNAL_SRFLUX'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIURNAL_SRFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPONGE
      MPI_master_only write(stdout,'(10x,A)') 'SPONGE'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'CLIMATOLOGY'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ZCLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'ZCLIMATOLOGY'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ZCLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2CLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'M2CLIMATOLOGY'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2CLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M3CLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'M3CLIMATOLOGY'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3CLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TCLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'TCLIMATOLOGY'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TCLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ZNUDGING
      MPI_master_only write(stdout,'(10x,A)') 'ZNUDGING'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ZNUDGING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2NUDGING
      MPI_master_only write(stdout,'(10x,A)') 'M2NUDGING'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2NUDGING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M3NUDGING
      MPI_master_only write(stdout,'(10x,A)') 'M3NUDGING'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3NUDGING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TNUDGING
      MPI_master_only write(stdout,'(10x,A)') 'TNUDGING'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TNUDGING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ROBUST_DIAG
      MPI_master_only write(stdout,'(10x,A)') 'ROBUST_DIAG'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ROBUST_DIAG'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'FRC_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Z_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'Z_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'M2_FRC_BRY'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M3_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'M3_FRC_BRY'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef T_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'T_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='T_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BSFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BSFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BSFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BTFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BTFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BTFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BMFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BMFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BMFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_INITIAL
      MPI_master_only write(stdout,'(10x,A)') 'ANA_INITIAL'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_INITIAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_SRFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_SRFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SRFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BRY
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ATM_PRESS
      MPI_master_only write(stdout,'(10x,A)') 'ATM_PRESS'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ATM_PRESS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PSOURCE
      MPI_master_only write(stdout,'(10x,A)') 'PSOURCE'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PSOURCE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_PSOURCE
      MPI_master_only write(stdout,'(10x,A)') 'ANA_PSOURCE'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_PSOURCE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M2FLATHER
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M2FLATHER'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2FLATHER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M2SPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M2SPECIFIED'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2SPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M2CHARACT
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M2CHARACT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2CHARACT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M2ORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M2ORLANSKI'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2ORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_VOLCONS
      MPI_master_only write(stdout,'(10x,A)') 'OBC_VOLCONS'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_VOLCONS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M3ORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M3ORLANSKI'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M3ORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_TORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'OBC_TORLANSKI'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_TORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M3SPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M3SPECIFIED'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M3SPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_TSPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'OBC_TSPECIFIED'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_TSPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AVERAGES
      MPI_master_only write(stdout,'(10x,A)') 'AVERAGES'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AVERAGES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AVERAGES_K
      MPI_master_only write(stdout,'(10x,A)') 'AVERAGES_K'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AVERAGES_K'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIAGNOSTICS_TS
      MPI_master_only write(stdout,'(10x,A)') 'DIAGNOSTICS_TS'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIAGNOSTICS_TS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIAGNOSTICS_UV
      MPI_master_only write(stdout,'(10x,A)') 'DIAGNOSTICS_UV'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIAGNOSTICS_UV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIAGNOSTICS_TS_ADV
      MPI_master_only write(stdout,'(10x,A)') 'DIAGNOSTICS_TS_ADV'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIAGNOSTICS_TS_ADV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIAGNOSTICS_TS_MLD
      MPI_master_only write(stdout,'(10x,A)') 'DIAGNOSTICS_TS_MLD'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIAGNOSTICS_TS_MLD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PISCES
      MPI_master_only write(stdout,'(10x,A)') 'PISCES'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PISCES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BIO_NChlPZD
      MPI_master_only write(stdout,'(10x,A)') 'BIO_NChlPZD'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BIO_NChlPZD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BIO_N2P2Z2D2
      MPI_master_only write(stdout,'(10x,A)') 'BIO_N2P2Z2D2'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BIO_N2P2Z2D2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BIO_N2ChlPZD2
      MPI_master_only write(stdout,'(10x,A)') 'BIO_N2ChlPZD2'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BIO_N2ChlPZD2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef key_trc_pisces
      MPI_master_only write(stdout,'(10x,A)') 'key_trc_pisces'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='key_trc_pisces'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef key_passivetrc
      MPI_master_only write(stdout,'(10x,A)') 'key_passivetrc'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='key_passivetrc'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIAGNOSTICS_BIO
      MPI_master_only write(stdout,'(10x,A)') 'DIAGNOSTICS_BIO'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIAGNOSTICS_BIO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef key_trc_diaadd
      MPI_master_only write(stdout,'(10x,A)') 'key_trc_diaadd'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='key_trc_diaadd'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef key_trc_dia3d
      MPI_master_only write(stdout,'(10x,A)') 'key_trc_dia3d'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='key_trc_dia3d'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OXYGEN
      MPI_master_only write(stdout,'(10x,A)') 'OXYGEN'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OXYGEN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VAR_CHL_C
      MPI_master_only write(stdout,'(10x,A)') 'VAR_CHL_C'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VAR_CHL_C'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef FLOATS_GLOBAL_ATTRIBUTES
      MPI_master_only write(stdout,'(10x,A)') 'FLOATS_GLOBAL_ATTRIBUTES'
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FLOATS_GLOBAL_ATTRIBUTES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef IBM
      MPI_master_only write(stdout,'(10x,A)') 'IBM'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='IBM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RANDOM_WALK
      MPI_master_only write(stdout,'(10x,A)') 'RANDOM_WALK'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RANDOM_WALK'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIEL_MIGRATION
      MPI_master_only write(stdout,'(10x,A)') 'DIEL_MIGRATION'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIEL_MIGRATION'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RANDOM_VERTICAL
      MPI_master_only write(stdout,'(10x,A)') 'RANDOM_VERTICAL'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RANDOM_VERTICAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RANDOM_HORIZONTAL
      MPI_master_only write(stdout,'(10x,A)') 'RANDOM_HORIZONTAL'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RANDOM_HORIZONTAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ALL_SIGMA
      MPI_master_only write(stdout,'(10x,A)') 'ALL_SIGMA'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ALL_SIGMA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_SEDIMENT
      MPI_master_only write(stdout,'(10x,A)') 'ANA_SEDIMENT'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SEDIMENT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BED_ARMOR
      MPI_master_only write(stdout,'(10x,A)') 'BED_ARMOR'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BED_ARMOR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_SPFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_SPFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SPFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BPFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BPFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BPFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LINEAR_CONTINUATION
      MPI_master_only write(stdout,'(10x,A)') 'LINEAR_CONTINUATION'
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LINEAR_CONTINUATION'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NEUMANN
      MPI_master_only write(stdout,'(10x,A)') 'NEUMANN'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NEUMANN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_WWAVE
      MPI_master_only write(stdout,'(10x,A)') 'ANA_WWAVE'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_WWAVE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BSEDIM
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BSEDIM'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BSEDIM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Z0_BL
      MPI_master_only write(stdout,'(10x,A)') 'Z0_BL'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z0_BL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Z0_RIP
      MPI_master_only write(stdout,'(10x,A)') 'Z0_RIP'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z0_RIP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Z0_BIO
      MPI_master_only write(stdout,'(10x,A)') 'Z0_BIO'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z0_BIO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ALLOW_SINGLE_BLOCK_MODE
      MPI_master_only write(stdout,'(10x,A)') 'ALLOW_SINGLE_BLOCK_MODE'
      ie=is +22
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ALLOW_SINGLE_BLOCK_MODE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SINGLE
      MPI_master_only write(stdout,'(10x,A)') 'SINGLE'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SINGLE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VAR_RHO_2D
      MPI_master_only write(stdout,'(10x,A)') 'VAR_RHO_2D'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VAR_RHO_2D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2FILTER_COSINE
      MPI_master_only write(stdout,'(10x,A)') 'M2FILTER_COSINE'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_COSINE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2FILTER_FLAT
      MPI_master_only write(stdout,'(10x,A)') 'M2FILTER_FLAT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_FLAT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2FILTER_POWER
      MPI_master_only write(stdout,'(10x,A)') 'M2FILTER_POWER'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_POWER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_SPLIT_UP3
      MPI_master_only write(stdout,'(10x,A)') 'UV_SPLIT_UP3'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_SPLIT_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef HADV_C4_UV
      MPI_master_only write(stdout,'(10x,A)') 'HADV_C4_UV'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HADV_C4_UV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_VIS4
      MPI_master_only write(stdout,'(10x,A)') 'UV_VIS4'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VIS4'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MIX_S_UV
      MPI_master_only write(stdout,'(10x,A)') 'MIX_S_UV'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MIX_S_UV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VIS_COEF_3D
      MPI_master_only write(stdout,'(10x,A)') 'VIS_COEF_3D'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VIS_COEF_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CLIMAT_UV_MIXH
      MPI_master_only write(stdout,'(10x,A)') 'CLIMAT_UV_MIXH'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CLIMAT_UV_MIXH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VADV_SPLINES_UV
      MPI_master_only write(stdout,'(10x,A)') 'VADV_SPLINES_UV'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_SPLINES_UV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VADV_C2_UV
      MPI_master_only write(stdout,'(10x,A)') 'VADV_C2_UV'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_C2_UV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef HADV_UPSTREAM_TS
      MPI_master_only write(stdout,'(10x,A)') 'HADV_UPSTREAM_TS'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HADV_UPSTREAM_TS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_DIF4
      MPI_master_only write(stdout,'(10x,A)') 'TS_DIF4'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_DIF4'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIF_COEF_3D
      MPI_master_only write(stdout,'(10x,A)') 'DIF_COEF_3D'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIF_COEF_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIF_SMAGO
      MPI_master_only write(stdout,'(10x,A)') 'DIF_SMAGO'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIF_SMAGO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef HADV_AKIMA_TS
      MPI_master_only write(stdout,'(10x,A)') 'HADV_AKIMA_TS'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HADV_AKIMA_TS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CLIMAT_TS_MIXH
      MPI_master_only write(stdout,'(10x,A)') 'CLIMAT_TS_MIXH'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CLIMAT_TS_MIXH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VADV_SPLINES_TS
      MPI_master_only write(stdout,'(10x,A)') 'VADV_SPLINES_TS'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_SPLINES_TS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VADV_AKIMA_TS
      MPI_master_only write(stdout,'(10x,A)') 'VADV_AKIMA_TS'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_AKIMA_TS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VADV_C2_TS
      MPI_master_only write(stdout,'(10x,A)') 'VADV_C2_TS'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_C2_TS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPONGE_DIF2
      MPI_master_only write(stdout,'(10x,A)') 'SPONGE_DIF2'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE_DIF2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPONGE_VIS2
      MPI_master_only write(stdout,'(10x,A)') 'SPONGE_VIS2'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE_VIS2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CONST_TRACERS
      MPI_master_only write(stdout,'(10x,A)') 'CONST_TRACERS'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CONST_TRACERS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_SM_UPDATE
      MPI_master_only write(stdout,'(10x,A)') 'BULK_SM_UPDATE'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_SM_UPDATE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DBLEPREC
      MPI_master_only write(stdout,'(10x,A)') 'DBLEPREC'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DBLEPREC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Linux
      MPI_master_only write(stdout,'(10x,A)') 'Linux'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Linux'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PGI
      MPI_master_only write(stdout,'(10x,A)') 'PGI'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PGI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef __IFC
      MPI_master_only write(stdout,'(10x,A)') '__IFC'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='__IFC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef QUAD
      MPI_master_only write(stdout,'(10x,A)') 'QUAD'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='QUAD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef QuadZero
      MPI_master_only write(stdout,'(10x,A)') 'QuadZero'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='QuadZero'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GLOBAL_2D_ARRAY
      MPI_master_only write(stdout,'(10x,A)') 'GLOBAL_2D_ARRAY'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLOBAL_2D_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GLOBAL_1D_ARRAYXI
      MPI_master_only write(stdout,'(10x,A)') 'GLOBAL_1D_ARRAYXI'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLOBAL_1D_ARRAYXI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GLOBAL_1D_ARRAYETA
      MPI_master_only write(stdout,'(10x,A)') 'GLOBAL_1D_ARRAYETA'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLOBAL_1D_ARRAYETA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef START_2D_ARRAY
      MPI_master_only write(stdout,'(10x,A)') 'START_2D_ARRAY'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='START_2D_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef START_1D_ARRAYXI
      MPI_master_only write(stdout,'(10x,A)') 'START_1D_ARRAYXI'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='START_1D_ARRAYXI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef START_1D_ARRAYETA
      MPI_master_only write(stdout,'(10x,A)') 'START_1D_ARRAYETA'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='START_1D_ARRAYETA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef EW_PERIODIC
      MPI_master_only write(stdout,'(10x,A)') 'EW_PERIODIC'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='EW_PERIODIC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NS_PERIODIC
      MPI_master_only write(stdout,'(10x,A)') 'NS_PERIODIC'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NS_PERIODIC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PRIVATE_1D_SCRATCH_ARRAY
      MPI_master_only write(stdout,'(10x,A)') 'PRIVATE_1D_SCRATCH_ARRAY'
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PRIVATE_1D_SCRATCH_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PRIVATE_2D_SCRATCH_ARRAY
      MPI_master_only write(stdout,'(10x,A)') 'PRIVATE_2D_SCRATCH_ARRAY'
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PRIVATE_2D_SCRATCH_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PRIVATE_1DXI_SCRATCH_ARRAY
      MPI_master_only write(stdout,'(10x,A)') 'PRIVATE_1DXI_SCRATCH_ARRAY'
      ie=is +25
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PRIVATE_1DXI_SCRATCH_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PRIVATE_1DETA_SCRATCH_ARRAY
      MPI_master_only write(stdout,'(10x,A)') 'PRIVATE_1DETA_SCRATCH_ARRAY'
      ie=is +26
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PRIVATE_1DETA_SCRATCH_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WESTERN_EDGE
      MPI_master_only write(stdout,'(10x,A)') 'WESTERN_EDGE'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WESTERN_EDGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef EASTERN_EDGE
      MPI_master_only write(stdout,'(10x,A)') 'EASTERN_EDGE'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='EASTERN_EDGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SOUTHERN_EDGE
      MPI_master_only write(stdout,'(10x,A)') 'SOUTHERN_EDGE'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SOUTHERN_EDGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NORTHERN_EDGE
      MPI_master_only write(stdout,'(10x,A)') 'NORTHERN_EDGE'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NORTHERN_EDGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MYID
      MPI_master_only write(stdout,'(10x,A)') 'MYID'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MYID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MPI_master_only
      MPI_master_only write(stdout,'(10x,A)') 'MPI_master_only'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI_master_only'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef mynode
      MPI_master_only write(stdout,'(10x,A)') 'mynode'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='mynode'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ZEROTH_TILE
      MPI_master_only write(stdout,'(10x,A)') 'ZEROTH_TILE'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ZEROTH_TILE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SINGLE_TILE_MODE
      MPI_master_only write(stdout,'(10x,A)') 'SINGLE_TILE_MODE'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SINGLE_TILE_MODE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LF_AM_STEP
      MPI_master_only write(stdout,'(10x,A)') 'LF_AM_STEP'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LF_AM_STEP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AB_AM_STEP
      MPI_master_only write(stdout,'(10x,A)') 'AB_AM_STEP'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AB_AM_STEP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef FIRST_TIME_STEP
      MPI_master_only write(stdout,'(10x,A)') 'FIRST_TIME_STEP'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FIRST_TIME_STEP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef FIRST_2D_STEP
      MPI_master_only write(stdout,'(10x,A)') 'FIRST_2D_STEP'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FIRST_2D_STEP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NOT_LAST_2D_STEP
      MPI_master_only write(stdout,'(10x,A)') 'NOT_LAST_2D_STEP'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NOT_LAST_2D_STEP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef sgi
      MPI_master_only write(stdout,'(10x,A)') 'sgi'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='sgi'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SGI
      MPI_master_only write(stdout,'(10x,A)') 'SGI'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SGI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CVECTOR
      MPI_master_only write(stdout,'(10x,A)') 'CVECTOR'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CVECTOR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CSDOACROSS
      MPI_master_only write(stdout,'(10x,A)') 'CSDOACROSS'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CSDOACROSS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CAND
      MPI_master_only write(stdout,'(10x,A)') 'CAND'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CAND'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ENTER_CRITICAL_REGION
      MPI_master_only write(stdout,'(10x,A)') 'ENTER_CRITICAL_REGION'
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ENTER_CRITICAL_REGION'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef EXIT_CRITICAL_REGION
      MPI_master_only write(stdout,'(10x,A)') 'EXIT_CRITICAL_REGION'
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='EXIT_CRITICAL_REGION'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CSDISTRIBUTE_RESHAPE
      MPI_master_only write(stdout,'(10x,A)') 'CSDISTRIBUTE_RESHAPE'
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CSDISTRIBUTE_RESHAPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BLOCK_PATTERN
      MPI_master_only write(stdout,'(10x,A)') 'BLOCK_PATTERN'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BLOCK_PATTERN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BLOCK_CLAUSE
      MPI_master_only write(stdout,'(10x,A)') 'BLOCK_CLAUSE'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BLOCK_CLAUSE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef cray
      MPI_master_only write(stdout,'(10x,A)') 'cray'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='cray'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CRAY
      MPI_master_only write(stdout,'(10x,A)') 'CRAY'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SHARE
      MPI_master_only write(stdout,'(10x,A)') 'SHARE'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SHARE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LOCAL
      MPI_master_only write(stdout,'(10x,A)') 'LOCAL'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LOCAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PUT_GRID_INTO_RESTART
      MPI_master_only write(stdout,'(10x,A)') 'PUT_GRID_INTO_RESTART'
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PUT_GRID_INTO_RESTART'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PUT_GRID_INTO_HISTORY
      MPI_master_only write(stdout,'(10x,A)') 'PUT_GRID_INTO_HISTORY'
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PUT_GRID_INTO_HISTORY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PUT_GRID_INTO_AVERAGES
      MPI_master_only write(stdout,'(10x,A)') 'PUT_GRID_INTO_AVERAGES'
      ie=is +21
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PUT_GRID_INTO_AVERAGES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NF_FTYPE
      MPI_master_only write(stdout,'(10x,A)') 'NF_FTYPE'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NF_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef nf_get_att_FTYPE
      MPI_master_only write(stdout,'(10x,A)') 'nf_get_att_FTYPE'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_get_att_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef nf_put_att_FTYPE
      MPI_master_only write(stdout,'(10x,A)') 'nf_put_att_FTYPE'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_put_att_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef nf_get_var1_FTYPE
      MPI_master_only write(stdout,'(10x,A)') 'nf_get_var1_FTYPE'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_get_var1_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef nf_put_var1_FTYPE
      MPI_master_only write(stdout,'(10x,A)') 'nf_put_var1_FTYPE'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_put_var1_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef nf_get_vara_FTYPE
      MPI_master_only write(stdout,'(10x,A)') 'nf_get_vara_FTYPE'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_get_vara_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef nf_put_vara_FTYPE
      MPI_master_only write(stdout,'(10x,A)') 'nf_put_vara_FTYPE'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_put_vara_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OUT_DOUBLE
      MPI_master_only write(stdout,'(10x,A)') 'OUT_DOUBLE'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OUT_DOUBLE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NF_FOUT
      MPI_master_only write(stdout,'(10x,A)') 'NF_FOUT'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NF_FOUT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef fast_indx_out
      MPI_master_only write(stdout,'(10x,A)') 'fast_indx_out'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='fast_indx_out'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_UPDATE_MIX
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_UPDATE_MIX'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_UPDATE_MIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_UPDATE_DECAL
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_UPDATE_DECAL'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_UPDATE_DECAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_SPONGE
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_SPONGE'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_SPONGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_EAST
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_EAST'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_EAST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_WEST
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_WEST'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_WEST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_NORTH
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_NORTH'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_NORTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_SOUTH
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_SOUTH'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_SOUTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_FLUX_BC
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_FLUX_BC'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_FLUX_BC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_M2SPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_M2SPECIFIED'
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_M2SPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_M3SPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_M3SPECIFIED'
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_M3SPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_TSPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_TSPECIFIED'
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_TSPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_M3ORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_M3ORLANSKI'
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_M3ORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_TORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_TORLANSKI'
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_TORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SMAGORINSKY
      MPI_master_only write(stdout,'(10x,A)') 'SMAGORINSKY'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SMAGORINSKY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
      MPI_master_only write(stdout,'(/)')
      if (iexample.eq.0) then
        MPI_master_only write(stdout,'(1x,A)')
     & 'ERROR in "cppdefs.h": no configuration is specified.'
        ierr=ierr+1
      elseif (iexample.gt.1) then
        MPI_master_only write(stdout,'(1x,A)')
     & 'ERROR: more than one configuration in "cppdefs.h".'
        ierr=ierr+1
      endif
      return
  99  MPI_master_only write(stdout,'(/1x,A,A/14x,A)')
     &  'CHECKDEFS -- ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
