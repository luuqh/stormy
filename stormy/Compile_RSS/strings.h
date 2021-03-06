! $Id: strings.h 697 2011-04-11 12:35:17Z gcambon $
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
! A long character string to hold activated cpp-switches.
! Basically it is used to keep track of cpp-switches by placing
! them together and writing into history file.
!                                    !
      integer max_opt_size           ! NOTE: Parameter max_opt_size
      parameter (max_opt_size=2048)  ! must be equal to the length
      character*2048 Coptions,srcs   ! of character string. 
      common /strings/ Coptions,srcs !
