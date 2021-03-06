! $Id: uv3dmix.F 697 2011-04-11 12:35:17Z gcambon $
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
#include "cppdefs.h"
#if defined SOLVE3D && defined UV_VIS2 && defined MIX_S_UV
# include "uv3dmix_S.F"
#elif defined SOLVE3D && defined UV_VIS2 && defined MIX_GP_UV
# include "uv3dmix_GP.F"
#elif defined SOLVE3D && defined UV_VIS4 && defined MIX_S_UV
# include "uv3dmix4_S.F"
#elif defined SOLVE3D && defined UV_VIS4 && defined MIX_GP_UV
# include "uv3dmix4_GP.F"
#else
      subroutine uv3dmix_empty
      end
#endif

