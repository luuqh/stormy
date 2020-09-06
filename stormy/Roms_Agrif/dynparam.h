! $Id: dynparam.h 697 2011-04-11 12:35:17Z gcambon $
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
! Dimensions of Physical Grid and array dimensions: 
! =========== == ======== ==== === ===== ============
! LLm,MMm  Number of the internal points of the PHYSICAL grid.
!          in the XI- and ETA-directions [physical side boundary
!          points and peroodic ghost points (if any) are excluded].
!
! Lm,Mm    Number of the internal points [see above] of array
!          covering a Message Passing subdomain. In the case when
!          no Message Passing partitioning is used, these two are
!          the same as LLm,MMm. 
!
      LLm=LLm0;  MMm=MMm0      ! LLm0,MMm0 defined in param.h
