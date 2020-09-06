! $Id: get_slp.F 697 2011-04-11 12:35:17Z Bijoy $
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
#ifdef ATM_PRESS

                              ! Read point or grided shortwave
      subroutine get_slp   ! radiation flux at the appropriate
                              ! time from forcing NetCDF file.
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "ncscrum.h"
#include "netcdf.inc"
      real cff
      integer i,ierr, lstr,lvar, lenstr, nf_fread, advance_cycle
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT
      if (itslp.eq.0 .or. iic.eq.0) then
        lstr=lenstr(frcname)
        if (may_day_flag.ne.0) return    !-->  EXIT
!

! If not opened yet, open forcing NetCDF file for reading.
! Find and save IDs for relevant variables, determine whether
! sea level pressure is a field or scalar value. 
!
        if (ncidfrc.eq.-1) then
          ierr=nf_open(frcname(1:lstr), nf_nowrite, ncidfrc)
          if (ierr .ne. nf_noerr) goto 4                !--> ERROR 
        endif

        ierr=nf_inq_varid (ncidfrc, 'sms_time', slp_tid)
        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'slp_time', frcname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!        write(stdout,3) 'slp time ==,', slp_time 

        lvar=lenstr(vname(1,indxSLP))
        ierr=nf_inq_varid (ncidfrc,vname(1,indxSLP)(1:lvar),slp_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidfrc, slp_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lslpgrd=1
            else
              lslpgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxSLP)(1:lvar),frcname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
! Determine whether there is cycling to reuse the input data
! and find cycling period "slp_cycle", set initial cycling
! index "slp_ncycle" and record index "slp_rec".
! Set initial value for time index "itslp" and both time record
! bounds to large negative artificial values, so that it will
! trigger the logic in reading part below.
!
        call set_cycle (ncidfrc,   slp_tid,    ntslp,
     &                  slp_cycle, slp_ncycle, slp_rec)
        if (may_day_flag.ne.0) return    !-->  EXIT
        itslp=2
        sms_time(1)=-1.E+20
        sms_time(2)=-1.E+20
        slp_scale=1.
      endif                   !<-- itslp.eq.0 .or. iic.eq.0

!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
  1   i=3-itslp
       cff=time+0.5*dt
        if (sms_time(i).le.cff .and. cff.lt.sms_time(itslp)) return
        ierr=advance_cycle (slp_cycle, ntslp, slp_ncycle, slp_rec)
        if (ierr.ne.0) goto 7                           !--> ERROR

        ierr=nf_get_var1_FTYPE(ncidfrc, slp_tid, slp_rec, cff)
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'slp_time', slp_rec 
          goto 99                                         !--> ERROR
        endif
        sms_time(i)=cff*day2sec+slp_cycle*slp_ncycle
        if (sms_time(itslp).eq.-1.E+20) sms_time(itslp)=sms_time(i)



        if (lslpgrd.eq.1) then
          ierr=nf_fread (slpg(START_2D_ARRAY,i), ncidfrc, slp_id,
     &                                             slp_rec, r2dvar)
        else
          ierr=nf_get_var1_FTYPE (ncidfrc,slp_id,slp_rec,slpp(i))
        endif

        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'slp', slp_rec
          goto 99                                         !--> ERROR
        endif      
       itslp=i
       MPI_master_only write(stdout,'(6x,A,1x,A,1x,g12.4,1x,I4)') 
     &           'GET_SLP --',
     &           'Read surface pressure for time =', cff
#ifdef MPI
     &                                                      , mynode
#endif
      if (ntslp.gt.1) goto 1
      if (ntslp.eq.1) return 
!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_SLP - unable to find forcing variable: ',a,
     &                          /,15x,'in forcing NetCDF file: ',a)
  4   write(stdout,5) frcname(1:lstr)
  5   format(/,' GET_SLP - unable to open forcing NetCDF file:',
     &                                                        1x,a)
      goto 99

  6   format(/,' GET_SLP - error while reading variable: ',a,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) slp_rec, ntslp, frcname(1:lstr), tdays,
     &                                sms_time(itslp)*sec2day
  8   format(/,' GET_SLP - ERROR: requested time record ',I4,
     &       1x,'exeeds the last available', /,14x,'record ',I4,
     &       1x,'in forcing NetCDF file: ',a, /,14x,'TDAYS = ',
     &       g12.4,2x,'last available SLP_TIME = ',g12.4)

  99  may_day_flag=2
      return
      end

      subroutine set_slp_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set-up surface pressure for current tile.
!--------------------------------------------------------------------
!
      implicit none
# include "param.h"
# include "forces.h"
# include "scalars.h"
# include "grid.h"
      integer Istr,Iend,Jstr,Jend, i,j, it1,it2
      real cff, cff1, cff2
      real val1
!
# include "compute_extended_bounds.h"
!
! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
! After that either load time-invariant data, or interpolate in time
! or complain about error and signal to quit, if interpolation is
! needed, but not possible.
!
!
      it1=3-itslp
      it2=itslp
      cff=time+0.5*dt
      cff1=sms_time(it2)-cff
      cff2=cff-sms_time(it1)
      if (ZEROTH_TILE.and. cff1.lt.dt) synchro_flag=.TRUE.
!
! Load time-invariant surface pressure.
! Time interpolation is not performed in this case.
!
      if (slp_cycle.lt.0.) then
        if (FIRST_TIME_STEP) then
          if (lslpgrd.eq.1) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                slp(i,j)=slp_scale*slpg(i,j,itslp)
              enddo
            enddo
          else
            val1=slp_scale*slpp(itslp)
            do j=JstrR,JendR
              do i=IstrR,IendR
                slp(i,j)=val1
              enddo
            enddo
          endif
        endif
!
! Time-interpolate shortwave radiation flux from grided or point
! data. Check that for the next time step [when time=time+dt]
! time+dt is still between slp_time(it1) and slp_time(it2);
! and if not, set synchro_flag to signal that the new forcing data
! should be read from the netCDF input file (master thread only).
!
      elseif (cff1.ge.0. .and. cff2.ge.0.) then
        cff=slp_scale/(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff
        if (lslpgrd.eq.1) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              slp(i,j)=cff1*slpg(i,j,it1)+cff2*slpg(i,j,it2)
            enddo
          enddo
        else
          val1=cff1*slpp(it1)+cff2*slpp(it2)
            do j=JstrR,JendR
              do i=IstrR,IendR
              slp(i,j)=val1
            enddo
          enddo
        endif
!
! Unable to set-up shortwave radiation flux:
! Complain about the error and signal to quit (ONE THREAD ONLY).
!
      else
        if (ZEROTH_TILE) then
          write(stdout,1) 'slp_time', tdays, sms_time(it2)*sec2day  
  1       format(/,' SET_SLP_TILE - current model time exceeds',
     &             ' ending value for variable: ',a,/,14x,'TDAYS = ',
     &                                      g12.4,2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif

      return
      end
#else
      subroutine get_slp_empty
      return
      end
#endif /* ATM_PRESS */
