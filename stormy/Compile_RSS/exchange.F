! $Id: exchange.F 697 2011-04-11 12:35:17Z gcambon $
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
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI

# define exchange_2d_tile exchange_r2d_tile
# define JSTART JstrR
# define ISTART IstrR
# include "exchange_2d_tile.h"
# undef ISTART 
# undef JSTART 
# undef exchange_2d_tile

# define exchange_2d_tile exchange_u2d_tile
# define JSTART JstrR
# define ISTART Istr
# include "exchange_2d_tile.h"
# undef ISTART 
# undef JSTART 
# undef exchange_2d_tile

# define exchange_2d_tile exchange_v2d_tile
# define JSTART Jstr
# define ISTART IstrR
# include "exchange_2d_tile.h"
# undef ISTART 
# undef JSTART 
# undef exchange_2d_tile

# define exchange_2d_tile exchange_p2d_tile
# define JSTART Jstr
# define ISTART Istr
# include "exchange_2d_tile.h"
# undef ISTART 
# undef JSTART 
# undef exchange_2d_tile

# ifdef SOLVE3D
#  define KSTART 1
#  define exchange_3d_tile exchange_r3d_tile
#  define JSTART JstrR
#  define ISTART IstrR
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile

#  define exchange_3d_tile exchange_u3d_tile
#  define JSTART JstrR
#  define ISTART Istr
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile

#  define exchange_3d_tile exchange_v3d_tile
#  define JSTART Jstr
#  define ISTART IstrR
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile

#  define exchange_3d_tile exchange_p3d_tile
#  define JSTART Jstr
#  define ISTART Istr
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile
#  undef KSTART

#  define KSTART 0
#  define exchange_3d_tile exchange_w3d_tile
#  define JSTART JstrR
#  define ISTART IstrR
#  include "exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef exchange_3d_tile
#  undef KSTART
# endif
#else
      subroutine exchange_empty
      end
#endif




