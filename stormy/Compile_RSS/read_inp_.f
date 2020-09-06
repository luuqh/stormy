      subroutine read_inp (ierr)
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=1519, MMm0=1399, N=32)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi/ Lmmpi,Mmmpi,
     &                    iminmpi,imaxmpi,jminmpi,jmaxmpi
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=1, NP_ETA=8,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6, Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      integer*4 NWEIGHT
      parameter (NWEIGHT=137)
      real dt, dtfast, time, time_start, tdays
      integer*4 iic, kstp, krhs, knew, next_kstp
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time,time_start, tdays,
     &                       iic, kstp, krhs, knew, next_kstp,
     &                       PREDICTOR_2D_STEP
      real time_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zob
     &               , xl, el, visc2, visc4, gamma2
      real  x_sponge,   v_sponge
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , ndtfast,nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
      logical ldefhis
      common /scalars_main/
     &             time_avg,  rho0,      rdrg,    rdrg2
     &           , Zob,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , ndtfast,nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , ldefhis
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers/CPU_time,proc,trd_count
		logical EAST_INTER
		logical WEST_INTER
		logical NORTH_INTER
		logical SOUTH_INTER
  	 integer mynode, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE, p_NW,p_NE
	 common /comm_setup/ mynode,ii,jj,p_W,p_E,p_S,p_N,p_SW,p_SE
	 common /comm_setup/ p_NW,p_NE, EAST_INTER, WEST_INTER
	 common /comm_setup/ NORTH_INTER, SOUTH_INTER
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  day2sec=86400.D0,
     &           sec2day=1.D0/86400.D0, jul_off=2440000.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxSSH
      parameter (indxSSH=indxVb+1)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+1, indxSVSTR=indxSSH+2)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+21)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+22)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+23)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+24)
	  integer indxSLP
	   parameter (indxSLP = indxSUSTR+25)
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      parameter (xi_rho=LLm+2,  xi_u=xi_rho-1,
     &           eta_rho=MMm+2, eta_v=eta_rho-1)
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms
     &      , ntsrf,  ntssh,  ntsst, ntsss, ntuclm, ntww,
     &        ntbulk
     &      , ntslp
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTstep, rstZ,    rstUb,  rstVb
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
      integer*4 ncidavg, nrecavg,  nrpfavg
     &      , avgTime, avgTstep, avgZ, avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUwstr, avgVwstr
      logical wrthis(14)
     &      , wrtavg(14)
      common/incscrum/
     &        ncidfrc, ncidbulk,ncidclm, ntsms, ntsrf, ntssh, ntsst
     &      , ntuclm, ntsss, ntww, ntbulk
     &      , ntslp
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTstep, rstZ,    rstUb,  rstVb
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , ncidavg,  nrecavg,  nrpfavg
     &      , avgTime, avgTstep, avgZ,    avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUWstr, avgVWstr
     &      , wrthis
     &      , wrtavg
      character*80 date_str, title
      character*80 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &                                ,   avgname
     &                                ,   bry_file
      character*52  vname(4, 39)
      common /cncscrum/       date_str,   title
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &                                ,  avgname
     &                                ,   bry_file
     &                      ,  vname
      include 'mpif.h'
      integer*4 kwsize, testunit, input
      parameter (kwsize=32, testunit=40, input=15)
      character end_signal*3, keyword*32, fname*64
      parameter (end_signal='end')
      integer*4 ierr, iargc, is,ie, kwlen, lstr, lenstr
      fname='roms.in'
      if (mynode.eq.0 .and. iargc().GT.0) call getarg(1,fname)
      call MPI_Bcast(fname,64,MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
      wrthis(indxTime)=.false.
      wrtavg(indxTime)=.false.
      ierr=0
      call setup_kwds (ierr)
      open (input,file=fname,status='old',form='formatted',err=97)
   1  keyword='                                '
      read(input,'(A)',err=1,end=99) keyword
      if (ichar(keyword(1:1)).eq.33) goto 1
      is=1
   2  if (is.eq.kwsize) then
        goto 1
      elseif (keyword(is:is).eq.' ') then
        is=is+1
        goto 2
      endif
      ie=is
   3  if (keyword(ie:ie).eq.':') then
        keyword(ie:ie)=' '
        goto 4
      elseif (keyword(ie:ie).ne.' ' .and. ie.lt.kwsize) then
        ie=ie+1
        goto 3
      endif
      goto 1
   4  kwlen=ie-is
      if (is.gt.1) keyword(1:kwlen)=keyword(is:is+kwlen-1)
      if (keyword(1:kwlen).eq.'title') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) title
        lstr=lenstr(title)
        if (mynode.eq.0) write(stdout,'(/1x,A)') title(1:lstr)
      elseif (keyword(1:kwlen).eq.'time_stepping') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ntimes,dt,ndtfast, ninfo
        if (mynode.eq.0) write(stdout,
     & '(I10,2x,A,1x,A /F10.2,2x,A,2(/I10,2x,A,1x,A)/F10.4,2x,A)'
     &  ) ntimes,  'ntimes   Total number of timesteps for',
     &                                            '3D equations.',
     &    dt,      'dt       Timestep [sec] for 3D equations',
     &    ndtfast, 'ndtfast  Number of 2D timesteps within each',
     &                                                 '3D step.',
     &    ninfo,   'ninfo    Number of timesteps between',
     &                                     'runtime diagnostics.'
        dtfast=dt/float(ndtfast)
        if (NWEIGHT.lt.(2*ndtfast-1)) then
          write(stdout,'(a,i3)')
     &    ' Error - Number of 2D timesteps (2*ndtfast-1): ',
     &    2*ndtfast-1
          write(stdout,'(a,i3)')
     &    '           exceeds barotopic weight dimension: ',NWEIGHT
          goto 95
        endif
      elseif (keyword(1:kwlen).eq.'initial') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrrec
        if (nrrec.gt.0) then
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          ininame=fname(1:lstr)
          if (mynode.eq.0) write(stdout,'(1x,A,2x,A,4x,A,I3)')
     &     'Initial State File:', ininame(1:lstr), 'Record:',nrrec
        endif
      elseif (keyword(1:kwlen).eq.'grid') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        open(testunit,file=fname(1:lstr), status='old', err=97)
        close(testunit)
        grdname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,'(10x,A,2x,A)')
     &                   'Grid File:', grdname(1:lstr)
      elseif (keyword(1:kwlen).eq.'forcing') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        open (testunit, file=fname(1:lstr), status='old', err=97)
        close(testunit)
        frcname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,'(2x,A,2x,A)')
     &             'Forcing Data File:', frcname(1:lstr)
      elseif (keyword(1:kwlen).eq.'restart') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrst, nrpfrst
        read(input,'(A)',err=95)  fname
        lstr=lenstr(fname)
        rstname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,
     &             '(7x,A,2x,A,4x,A,I6,4x,A,I4)')
     &             'Restart File:', rstname(1:lstr),
     &             'nrst =', nrst, 'rec/file: ', nrpfrst
      elseif (keyword(1:kwlen).eq.'history') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefhis, nwrt, nrpfhis
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        hisname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,
     &             '(7x,A,2x,A,2x,A,1x,L1,2x,A,I4,2x,A,I3)')
     &       'History File:', hisname(1:lstr),  'Create new:',
     &       ldefhis, 'nwrt =', nwrt, 'rec/file =', nrpfhis
      elseif (keyword(1:kwlen).eq.'averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ntsavg, navg, nrpfavg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        avgname=fname(1:lstr)
        if (mynode.eq.0) write(stdout,
     &         '(2(I10,2x,A,1x,A/32x,A/),6x,A,2x,A,1x,A,I3)')
     &      ntsavg, 'ntsavg      Starting timestep for the',
     &         'accumulation of output', 'time-averaged data.',
     &      navg,   'navg        Number of timesteps between',
     &     'writing of time-averaged','data into averages file.',
     &     'Averages File:', avgname(1:lstr),
     &     'rec/file =', nrpfavg
      elseif (keyword(1:kwlen).eq.'primary_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrthis(indxZ),  wrthis(indxUb)
     &                                       ,  wrthis(indxVb)
        if ( wrthis(indxZ) .or. wrthis(indxUb) .or. wrthis(indxVb)
     &     ) wrthis(indxTime)=.true.
        if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &    'Fields to be saved in history file: (T/F)'
     &    , wrthis(indxZ),  'write zeta ', 'free-surface.'
     &    , wrthis(indxUb), 'write UBAR ', '2D U-momentum component.'
     &    , wrthis(indxVb), 'write VBAR ', '2D V-momentum component.'
      elseif (keyword(1:kwlen).eq.'primary_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxZ),  wrtavg(indxUb)
     &                                    ,  wrtavg(indxVb)
        if ( wrtavg(indxZ) .or. wrtavg(indxUb) .or. wrtavg(indxVb)
     &     ) wrtavg(indxTime)=.true.
        if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &  'Fields to be saved in averages file: (T/F)'
     &  , wrtavg(indxZ),  'write zeta ', 'free-surface.'
     &  , wrtavg(indxUb), 'write UBAR ', '2D U-momentum component.'
     &  , wrtavg(indxVb), 'write VBAR ', '2D V-momentum component.'
      elseif (keyword(1:kwlen).eq.'rho0') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rho0
        if (mynode.eq.0) write(stdout,'(F10.4,2x,A,1x,A)')
     &        rho0, 'rho0     Boussinesq approximation',
     &                           'mean density, kg/m3.'
      elseif (keyword(1:kwlen).eq.'lateral_visc') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) visc2, visc4
        if (mynode.eq.0) write(stdout,9) visc2
   9    format(1pe10.3,2x,'visc2    Horizontal Laplacian ',
     &       'mixing coefficient [m2/s]',/,32x,'for momentum.')
      elseif (keyword(1:kwlen).eq.'bottom_drag') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rdrg, rdrg2, Zob, Cdb_min, Cdb_max
        if (mynode.eq.0) write(stdout,'(5(1pe10.3,2x,A/))')
     &     rdrg, 'rdrg     Linear bottom drag coefficient (m/si).',
     &    rdrg2, 'rdrg2    Quadratic bottom drag coefficient.',
     &      Zob, 'Zob      Bottom roughness for logarithmic law (m).',
     &  Cdb_min, 'Cdb_min  Minimum bottom drag coefficient.',
     &  Cdb_max, 'Cdb_max  Maximum bottom drag coefficient.'
      elseif (keyword(1:kwlen).eq.'gamma2') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) gamma2
        if (mynode.eq.0) write(stdout,'(f10.2,2x,A,1x,A)')
     &     gamma2, 'gamma2   Slipperiness parameter:',
     &                     'free-slip +1, or no-slip -1.'
      elseif (keyword(1:kwlen).eq.'sponge') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) x_sponge, v_sponge
        if (mynode.eq.0) write(stdout,'(1pe10.2,2x,A,1x,A)')
     &     x_sponge,'x_sponge Thickness of sponge',
     &     'and/or nudging layer (m)'
        if (mynode.eq.0) write(stdout,'(f10.2,2x,A)')
     &     v_sponge,'v_sponge Viscosity in sponge layer (m2/s)'
      elseif (keyword(1:kwlen).eq.'nudg_cof') then
        call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) tauT_in,tauT_out,tauM_in,tauM_out
          tauT_in =1.D0/(tauT_in *86400.D0)
          tauT_out=1.D0/(tauT_out*86400.D0)
          tauM_in =1.D0/(tauM_in *86400.D0)
          tauM_out=1.D0/(tauM_out*86400.D0)
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A)')
     &        tauT_in,'tauT_in  Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A)')
     &       tauT_out,'tauT_out Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A)')
     &        tauM_in,'tauM_in  Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A/)')
     &       tauM_out,'tauM_out Nudging coefficients [sec^-1]'
      else
        if (mynode.eq.0) write(stdout,'(/3(1x,A)/)')
     &                  'WARNING: Unrecognized keyword:',
     &                   keyword(1:kwlen),' --> DISREGARDED.'
      endif
      if (keyword(1:kwlen) .eq. end_signal) goto 99
      goto 1
  95  write(stdout,'(/1x,4A/)') 'READ_INP ERROR while reading block',
     &                    ' with keyword ''', keyword(1:kwlen), '''.'
      ierr=ierr+1
      goto 99
  97  write(stdout,'(/1x,4A/)') 'READ_INP ERROR: Cannot find input ',
     &                                'file ''', fname(1:lstr), '''.'
      ierr=ierr+1
  99  close (input)
      if (ierr.eq.0) then
        call check_kwds (ierr)
        call check_srcs
        call check_switches1 (ierr)
        call check_switches2 (ierr)
      endif
      if (ierr.ne.0) then
        write(stdout,'(/1x,2A,I3,1x,A/)') 'READ_INP ERROR: ',
     & 'A total of', ierr, 'configuration errors discovered.'
        return
      endif
      return
      end
