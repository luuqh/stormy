! $Id: work.h 697 2011-04-11 12:35:17Z gcambon $
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
!
! This is "work.h": declaration of utility work array.
!
#ifdef SOLVE3D
      real work(GLOBAL_2D_ARRAY,0:N)
      common /work3d/ work
#endif

      real work2d(GLOBAL_2D_ARRAY)
      common /work2d/ work2d

      real work2d2(GLOBAL_2D_ARRAY)
      common /work2d2/ work2d2
