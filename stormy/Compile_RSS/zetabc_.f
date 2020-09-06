      subroutine zetabc_tile(Istr,Iend,Jstr,Jend)
      implicit none
      integer*4 Istr,Iend,Jstr,Jend, i,j
      real cff,cx,dft,dfx, eps
      parameter (eps=1.D-20)
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
      real zetabry_west(-1:Mm+2+padd_E),
     &    zetabry_west_dt(-1:Mm+2+padd_E,2)
      common /bry_zeta_west/ zetabry_west, zetabry_west_dt
      real ubarbry_west(-1:Mm+2+padd_E),
     &    ubarbry_west_dt(-1:Mm+2+padd_E,2)
     &    ,vbarbry_west(-1:Mm+2+padd_E),
     &    vbarbry_west_dt(-1:Mm+2+padd_E,2)
      common /bry_ubar_west/ ubarbry_west, ubarbry_west_dt,
     &                       vbarbry_west, vbarbry_west_dt
      real zetabry_east(-1:Mm+2+padd_E),
     &    zetabry_east_dt(-1:Mm+2+padd_E,2)
      common /bry_zeta_east/ zetabry_east, zetabry_east_dt
      real ubarbry_east(-1:Mm+2+padd_E),
     &    ubarbry_east_dt(-1:Mm+2+padd_E,2)
     &    ,vbarbry_east(-1:Mm+2+padd_E),
     &    vbarbry_east_dt(-1:Mm+2+padd_E,2)
      common /bry_ubar_east/ ubarbry_east, ubarbry_east_dt,
     &                       vbarbry_east, vbarbry_east_dt
      real zetabry_south(-1:Lm+2+padd_X),
     &    zetabry_south_dt(-1:Lm+2+padd_X,2)
      common /bry_zeta_south/ zetabry_south, zetabry_south_dt
      real ubarbry_south(-1:Lm+2+padd_X),
     &    ubarbry_south_dt(-1:Lm+2+padd_X,2)
     &    ,vbarbry_south(-1:Lm+2+padd_X),
     &    vbarbry_south_dt(-1:Lm+2+padd_X,2)
      common /bry_ubar_south/ ubarbry_south, ubarbry_south_dt,
     &                        vbarbry_south, vbarbry_south_dt
      real zetabry_north(-1:Lm+2+padd_X),
     &    zetabry_north_dt(-1:Lm+2+padd_X,2)
      common /bry_zeta_north/ zetabry_north, zetabry_north_dt
      real ubarbry_north(-1:Lm+2+padd_X),
     &    ubarbry_north_dt(-1:Lm+2+padd_X,2)
     &    ,vbarbry_north(-1:Lm+2+padd_X),
     &    vbarbry_north_dt(-1:Lm+2+padd_X,2)
      common /bry_ubar_north/ ubarbry_north, ubarbry_north_dt,
     &                        vbarbry_north, vbarbry_north_dt
      real h(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real hinv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real f(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real fomn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real angler(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_angler/angler
      real latr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_latr/latr /grid_lonr/lonr
      real pm(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pm_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pm_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
     &       /metrics_omr/om_r /metrics_on_r/on_r
     &       /metrics_omu/om_u /metrics_on_u/on_u
     &       /metrics_omv/om_v /metrics_on_v/on_v
     &       /metrics_omp/om_p /metrics_on_p/on_p
     &       /metrics_pnu/pn_u /metrics_pmv/pm_v
     &       /metrics_pmu/pm_u /metrics_pnv/pn_v
      real dmde(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real dndx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_dmde/dmde    /metrics_dndx/dndx
      real pmon_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmon_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmon_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real grdscl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
     &       /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
     &       /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
     &                              /metrics_grdscl/grdscl
      real rmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real umask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real vmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mask_r/rmask /mask_p/pmask
     &       /mask_u/umask /mask_v/vmask
      real zeta(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real ubar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real vbar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      common /ocean_zeta/zeta
     &       /ocean_ubar/ubar
     &       /ocean_vbar/vbar
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
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      integer*4 JstrV
      if (.not.WEST_INTER) then
        IstrR=Istr-1
        IstrU=Istr+1
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if (.not.EAST_INTER) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if (.not.SOUTH_INTER) then
        JstrR=Jstr-1
        JstrV=Jstr+1
      else
        JstrR=Jstr
        JstrV=Jstr
      endif
      if (.not.NORTH_INTER) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      if (.not.WEST_INTER) then
        do j=JstrV-1,Jend
          cx=dtfast*pm(1,j)*sqrt(g*h(1,j))
          zeta(0,j,knew)=(zeta(0,j,kstp)+cx*zeta(1,j,knew))/(1.D0+cx)
     &                      *rmask(0,j)
        enddo
      endif
      if (.not.EAST_INTER) then
        do j=JstrV-1,Jend
          cx=dtfast*pm(Iend,j)*sqrt(g*h(Iend,j))
          zeta(Iend+1,j,knew)=(zeta(Iend+1,j,kstp)
     &                              +cx*zeta(Iend,j,knew))/(1.D0+cx)
     &                      *rmask(Iend+1,j)
        enddo
      endif
      if (.not.SOUTH_INTER) then
        do i=IstrU-1,Iend
          cx=dtfast*pn(i,1)*sqrt(g*h(i,1))
          zeta(i,0,knew)=(zeta(i,0,kstp)
     &                              +cx*zeta(i,1,knew))/(1.D0+cx)
     &                      *rmask(i,0)
        enddo
      endif
      if (.not.NORTH_INTER) then
        do i=IstrU-1,Iend
          cx=dtfast*pn(i,Jend)*sqrt(g*h(i,Jend))
          zeta(i,Jend+1,knew)=(zeta(i,Jend+1,kstp)
     &                              +cx*zeta(i,Jend,knew))/(1.D0+cx)
     &                      *rmask(i,Jend+1)
        enddo
      endif
      if (.not.SOUTH_INTER .and. .not.WEST_INTER) then
        zeta(0,0,knew)=0.5D0*(zeta(1,0 ,knew)+zeta(0 ,1,knew))
     &                       *rmask(0,0)
      endif
      if (.not.SOUTH_INTER .and. .not.EAST_INTER) then
        zeta(Iend+1,0,knew)=0.5D0*(zeta(Iend+1,1 
     &                                         ,knew)+zeta(Iend,0,knew))
     &                       *rmask(Iend+1,0)
      endif
      if (.not.NORTH_INTER .and. .not.WEST_INTER) then
        zeta(0,Jend+1,knew)=0.5D0*(zeta(0,Jend,knew)+zeta(1 
     &                                                    ,Jend+1,knew))
     &                       *rmask(0,Jend+1)
      endif
      if (.not.NORTH_INTER .and. .not.EAST_INTER) then
        zeta(Iend+1,Jend+1,knew)=0.5D0*(zeta(Iend+1,Jend,knew)+
     &                            zeta(Iend,Jend+1,knew))
     &                       *rmask(Iend+1,Jend+1)
      endif
      return
      end
