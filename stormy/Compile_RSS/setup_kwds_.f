      subroutine setup_kwds (ierr)
      implicit none
      integer*4 ierr, is,ie
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
      integer*4 max_opt_size
      parameter (max_opt_size=2048)
      character*2048 Coptions,srcs
      common /strings/ Coptions,srcs
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
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      is=1
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='title'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='time_stepping'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='initial'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='grid'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='forcing'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='restart'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='history'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='averages'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +22
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='primary_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='primary_averages'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='rho0'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='lateral_visc'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='bottom_drag'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='gamma2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='sponge'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nudg_cof'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      return
  99  if (mynode.eq.0) write(stdout,'(/1x,A,A/14x,A)')
     &  'SETUP_KWDS ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
