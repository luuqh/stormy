      subroutine MPI_Setup (ierr)
      integer*4 ierr, nsize, ii_W, ii_E, jj_S, jj_N
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istrmpi,Iendmpi,Jstrmpi,Jendmpi, i_X,j_E
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
      call MPI_Comm_size (MPI_COMM_WORLD, nsize,  ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, mynode, ierr)
      if (nsize.ne.NNODES) then
        if (mynode.eq.0) write(stdout,'(/1x,A,I4,1x,A,I3,A/)')
     &     'ERROR in MPI_Setup: number of MPI-nodes should be',
     &                         NNODES, 'instead of', nsize, '.'
        ierr=99
        return
      endif
      ii=mod(mynode,NP_XI)
      jj=mynode/NP_XI
      if (NP_XI.eq.1) then
        WEST_INTER=.false.
        EAST_INTER=.false.
      else
        if (ii.eq.0) then
          WEST_INTER=.false.
        else
          WEST_INTER=.true.
        endif
        if (ii.eq.NP_XI-1) then
          EAST_INTER=.false.
        else
          EAST_INTER=.true.
        endif
      endif
      if (NP_ETA.eq.1) then
        SOUTH_INTER=.false.
        NORTH_INTER=.false.
      else
        if (jj.eq.0) then
          SOUTH_INTER=.false.
        else
          SOUTH_INTER=.true.
        endif
        if (jj.eq.NP_ETA-1) then
          NORTH_INTER=.false.
        else
          NORTH_INTER=.true.
        endif
      endif
      ii_W=mod(ii-1+NP_XI,NP_XI)
      ii_E=mod(ii+1       ,NP_XI)
      jj_S=mod(jj-1+NP_ETA,NP_ETA)
      jj_N=mod(jj+1       ,NP_ETA)
      p_W=ii_W +NP_XI*jj
      p_E=ii_E +NP_XI*jj
      p_S=ii   +NP_XI*jj_S
      p_N=ii   +NP_XI*jj_N
      p_NW=ii_W+NP_XI*jj_N
      p_SW=ii_W+NP_XI*jj_S
      p_NE=ii_E+NP_XI*jj_N
      p_SE=ii_E+NP_XI*jj_S
      chunk_size_X=(LLm+NP_XI-1)/NP_XI
      margin_X=(NP_XI*chunk_size_X-LLm)/2
      chunk_size_E=(MMm+NP_ETA-1)/NP_ETA
      margin_E=(NP_ETA*chunk_size_E-MMm)/2
      j_E=mynode/NP_XI
      i_X=mynode-j_E*NP_XI
      istrmpi=1+i_X*chunk_size_X-margin_X
      iendmpi=istrmpi+chunk_size_X-1
      istrmpi=max(istrmpi,1)
      iendmpi=min(iendmpi,LLm)
      jstrmpi=1+j_E*chunk_size_E-margin_E
      jendmpi=jstrmpi+chunk_size_E-1
      jstrmpi=max(jstrmpi,1)
      jendmpi=min(jendmpi,MMm)
      Lmmpi=iendmpi-istrmpi+1
      Mmmpi=jendmpi-jstrmpi+1
      iminmpi=istrmpi
      imaxmpi=iendmpi
      jminmpi=jstrmpi
      jmaxmpi=jendmpi
      return
      end
