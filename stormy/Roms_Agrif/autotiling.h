! $Id: autotiling.h 697 2011-04-11 12:35:17Z gcambon $
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
         integer nbsampling, lastiic
         common/dyndistribution/nbsampling, lastiic

#if !defined AGRIF
      integer, parameter :: MAX_NSUB_X=Lm/6
      integer, parameter :: MAX_NSUB_E=Mm/6
#else
      integer :: MAX_NSUB_X
      integer :: MAX_NSUB_E
      common/autotilingparams/MAX_NSUB_X,MAX_NSUB_E
#endif
      
      integer, parameter :: nbvalid = 3
      integer, parameter :: nbdistrib = 4
      
      real cpu_domain(0:2)
      real times(MAX_NSUB_X,MAX_NSUB_E,nbdistrib*nbvalid)
      integer nbinst(MAX_NSUB_X,MAX_NSUB_E)
      integer nsub_emin,nsub_emax, nsub_xmin,nsub_xmax
      integer nsub_decalx, nsub_decale
      integer curdistrib
      common/dyndistribution2/times,cpu_domain,
     &       nbinst,nsub_xmin, nsub_xmax,
     &  nsub_emin,nsub_emax,nsub_decalx,nsub_decale,
     &  curdistrib

      integer, parameter :: nbx=4
      integer, parameter :: nbe=4

