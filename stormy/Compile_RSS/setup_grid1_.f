      subroutine setup_grid1 (tile)
      implicit none
      integer*4 tile, trd
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
C$    integer*4 omp_get_thread_num
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lmmpi+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lmmpi)/2
      chunk_size_E=(Mmmpi+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mmmpi)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lmmpi)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mmmpi)
      call setup_grid1_tile (Istr,Iend,Jstr,Jend)
      return
      end
      subroutine setup_grid1_tile (Istr,Iend,Jstr,Jend)
      implicit none
      integer*4 Istr,Iend,Jstr,Jend, i,j
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
      integer*4 IstrR,IendR,JstrR,JendR
      if (Istr.eq.1) then
        if (WEST_INTER) then
          IstrR=Istr-2
        else
          IstrR=Istr-1
        endif
      else
        IstrR=Istr
      endif
      if (Iend.eq.Lmmpi) then
        if (EAST_INTER) then
          IendR=Iend+2
        else
          IendR=Iend+1
        endif
      else
        IendR=Iend
      endif
      if (Jstr.eq.1) then
        if (SOUTH_INTER) then
          JstrR=Jstr-2
        else
          JstrR=Jstr-1
        endif
      else
        JstrR=Jstr
      endif
      if (Jend.eq.Mmmpi) then
        if (NORTH_INTER) then
          JendR=Jend+2
        else
          JendR=Jend+1
        endif
      else
        JendR=Jend
      endif
      do j=JstrR,JendR
        do i=IstrR,IendR
          fomn(i,j)=f(i,j)/(pm(i,j)*pn(i,j))
        enddo
      enddo
      if (WEST_INTER) IstrR=Istr
      if (EAST_INTER) IendR=Iend
      if (SOUTH_INTER) JstrR=Jstr
      if (NORTH_INTER) JendR=Jend
      do j=JstrR,JendR
        do i=IstrR,IendR
          om_r(i,j)=1.D0/pm(i,j)
          on_r(i,j)=1.D0/pn(i,j)
          pnom_r(i,j)=pn(i,j)/pm(i,j)
          pmon_r(i,j)=pm(i,j)/pn(i,j)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          dndx(i,j)=0.5D0/pn(i+1,j)-0.5D0/pn(i-1,j)
          dmde(i,j)=0.5D0/pm(i,j+1)-0.5D0/pm(i,j-1)
        enddo
      enddo
      do j=JstrR,JendR
        do i=Istr,IendR
           pmon_u(i,j)=(pm(i,j)+pm(i-1,j))
     &                 /(pn(i,j)+pn(i-1,j))
           om_u(i,j)=2.D0/(pm(i,j)+pm(i-1,j))
           on_u(i,j)=2.D0/(pn(i,j)+pn(i-1,j))
           pn_u(i,j)=0.5D0*(pn(i,j)+pn(i-1,j))
           pm_u(i,j)=0.5D0*(pm(i,j)+pm(i-1,j))
           umask(i,j)=rmask(i,j)*rmask(i-1,j)
        enddo
      enddo
      do j=Jstr,JendR
        do i=IstrR,IendR
          pnom_v(i,j)=(pn(i,j)+pn(i,j-1))
     &                /(pm(i,j)+pm(i,j-1))
          om_v(i,j)=2.D0/(pm(i,j)+pm(i,j-1))
          on_v(i,j)=2.D0/(pn(i,j)+pn(i,j-1))
          pm_v(i,j)=0.5D0*(pm(i,j)+pm(i,j-1))
          pn_v(i,j)=0.5D0*(pn(i,j)+pn(i,j-1))
          vmask(i,j)=rmask(i,j)*rmask(i,j-1)
        enddo
      enddo
      do j=Jstr,JendR
        do i=Istr,IendR
          pnom_p(i,j)=(pn(i,j)+pn(i,j-1)+pn(i-1,j)+pn(i-1,j-1))
     &               /(pm(i,j)+pm(i,j-1)+pm(i-1,j)+pm(i-1,j-1))
          pmon_p(i,j)=(pm(i,j)+pm(i,j-1)+pm(i-1,j)+pm(i-1,j-1))
     &               /(pn(i,j)+pn(i,j-1)+pn(i-1,j)+pn(i-1,j-1))
          om_p(i,j)=4.D0/(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          on_p(i,j)=4.D0/(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
          pmask(i,j)=rmask(i,j)*rmask(i-1,j)*rmask(i,j-1)
     &                                      *rmask(i-1,j-1)
          if (gamma2.lt.0.D0) pmask(i,j)=2.D0-pmask(i,j)
        enddo
      enddo
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,   om_r)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,   on_r)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend, pnom_r)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend, pmon_r)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,   dndx)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,   dmde)
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend, pmon_u)
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,   om_u)
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,   on_u)
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,   pn_u)
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,   pm_u)
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend, pnom_v)
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,   om_v)
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,   on_v)
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,   pm_v)
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,   pn_v)
      call exchange_p2d_tile (Istr,Iend,Jstr,Jend, pnom_p)
      call exchange_p2d_tile (Istr,Iend,Jstr,Jend, pmon_p)
      call exchange_p2d_tile (Istr,Iend,Jstr,Jend,   om_p)
      call exchange_p2d_tile (Istr,Iend,Jstr,Jend,   on_p)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,  rmask)
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,  umask)
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,  vmask)
      call exchange_p2d_tile (Istr,Iend,Jstr,Jend,  pmask)
      return
      end