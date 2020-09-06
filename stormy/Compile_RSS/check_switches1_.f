      subroutine check_switches1 (ierr)
      implicit none
      integer*4 ierr, is,ie, iexample
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
      if (mynode.eq.0) write(stdout,'(/1x,A/)')
     &      'Activated C-preprocessing Options:'
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      iexample=0
      is=1
      iexample=iexample+1
      if (mynode.eq.0) write(stdout,'(10x,A)') 'REGIONAL'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='REGIONAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'STSU'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='STSU'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'MPI'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_EAST'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_EAST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_WEST'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_WEST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_NORTH'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_NORTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_SOUTH'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_SOUTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_COR'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_COR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_ADV'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_ADV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'CURVGRID'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CURVGRID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SPHERICAL'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPHERICAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'MASKING'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MASKING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_VIS2'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VIS2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SPONGE'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'FRC_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'Z_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'M2_FRC_BRY'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_BSFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BSFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_BTFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BTFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_INITIAL'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_INITIAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_SRFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SRFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ATM_PRESS'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ATM_PRESS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_M2CHARACT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2CHARACT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'AVERAGES'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AVERAGES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'M2FILTER_POWER'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_POWER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'VADV_SPLINES_UV'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_SPLINES_UV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'HADV_UPSTREAM_TS'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HADV_UPSTREAM_TS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'VADV_AKIMA_TS'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_AKIMA_TS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'DBLEPREC'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DBLEPREC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'Linux'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Linux'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'QUAD'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='QUAD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'QuadZero'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='QuadZero'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLOBAL_2D_ARRAY'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLOBAL_2D_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLOBAL_1D_ARRAYXI'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLOBAL_1D_ARRAYXI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLOBAL_1D_ARRAYETA'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLOBAL_1D_ARRAYETA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'START_2D_ARRAY'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='START_2D_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'START_1D_ARRAYXI'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='START_1D_ARRAYXI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'START_1D_ARRAYETA'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='START_1D_ARRAYETA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 
     &                                        'PRIVATE_1D_SCRATCH_ARRAY'
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PRIVATE_1D_SCRATCH_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 
     &                                        'PRIVATE_2D_SCRATCH_ARRAY'
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PRIVATE_2D_SCRATCH_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 
     &                                      'PRIVATE_1DXI_SCRATCH_ARRAY'
      ie=is +25
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PRIVATE_1DXI_SCRATCH_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 
     &                                     'PRIVATE_1DETA_SCRATCH_ARRAY'
      ie=is +26
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PRIVATE_1DETA_SCRATCH_ARRAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'WESTERN_EDGE'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WESTERN_EDGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'EASTERN_EDGE'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='EASTERN_EDGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SOUTHERN_EDGE'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SOUTHERN_EDGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NORTHERN_EDGE'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NORTHERN_EDGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'MYID'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MYID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'MPI_master_only'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI_master_only'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ZEROTH_TILE'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ZEROTH_TILE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SINGLE_TILE_MODE'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SINGLE_TILE_MODE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'LF_AM_STEP'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LF_AM_STEP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'FIRST_TIME_STEP'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FIRST_TIME_STEP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'FIRST_2D_STEP'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FIRST_2D_STEP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NOT_LAST_2D_STEP'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NOT_LAST_2D_STEP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'PUT_GRID_INTO_RESTART'
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PUT_GRID_INTO_RESTART'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'PUT_GRID_INTO_HISTORY'
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PUT_GRID_INTO_HISTORY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'PUT_GRID_INTO_AVERAGES'
      ie=is +21
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PUT_GRID_INTO_AVERAGES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NF_FTYPE'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NF_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'nf_get_att_FTYPE'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_get_att_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'nf_put_att_FTYPE'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_put_att_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'nf_get_var1_FTYPE'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_get_var1_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'nf_put_var1_FTYPE'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_put_var1_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'nf_get_vara_FTYPE'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_get_vara_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'nf_put_vara_FTYPE'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nf_put_vara_FTYPE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NF_FOUT'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NF_FOUT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'fast_indx_out'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='fast_indx_out'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'AGRIF_UPDATE_MIX'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_UPDATE_MIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'AGRIF_UPDATE_DECAL'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_UPDATE_DECAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'AGRIF_SPONGE'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_SPONGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(/)')
      if (iexample.eq.0) then
        if (mynode.eq.0) write(stdout,'(1x,A)')
     & 'ERROR in "cppdefs.h": no configuration is specified.'
        ierr=ierr+1
      elseif (iexample.gt.1) then
        if (mynode.eq.0) write(stdout,'(1x,A)')
     & 'ERROR: more than one configuration in "cppdefs.h".'
        ierr=ierr+1
      endif
      return
  99  if (mynode.eq.0) write(stdout,'(/1x,A,A/14x,A)')
     &  'CHECKDEFS -- ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
