      subroutine step2d (tile)
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
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,4,0:NPP-1)
      common /private_scratch/ A2d,A3d
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
      trd=0
C$    trd=omp_get_thread_num()
      call step2D_FB_tile (istr,iend,jstr,jend, A2d(1,1,trd),
     &                    A2d(1,2,trd), A2d(1, 3,trd), A2d(1, 4,trd),
     &                    A2d(1, 5,trd), A2d(1, 6,trd),A2d(1, 7,trd),
     &                    A2d(1, 8,trd), A2d(1, 9,trd), A2d(1,10,trd),
     &                    A2d(1,11,trd), A2d(1,12,trd), A2d(1,13,trd))
      return
      end
      subroutine step2D_FB_tile (istr,iend,jstr,jend, zeta_new,Dnew,
     &                           rubar,rvbar, urhs,vrhs, DUon,DVom,
     &                           Drhs, UFx,UFe,VFx,VFe)
      implicit none
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
        real fac
      integer*4 istr,iend,jstr,jend, i,j, kbak, kold,
     &        IMAX,JMAX
      real    VMAX,VMAXL
      real zeta_new(Istr-2:Iend+2,Jstr-2:Jend+2),  cff,
     &         Dnew(Istr-2:Iend+2,Jstr-2:Jend+2),  cff0,
     &        rubar(Istr-2:Iend+2,Jstr-2:Jend+2),  cff1,
     &        rvbar(Istr-2:Iend+2,Jstr-2:Jend+2),  cff2,
     &         urhs(Istr-2:Iend+2,Jstr-2:Jend+2),  cff3,
     &         vrhs(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         DUon(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         DVom(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         Drhs(Istr-2:Iend+2,Jstr-2:Jend+2),
     &          UFx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &          UFe(Istr-2:Iend+2,Jstr-2:Jend+2),  DUnew,
     &          VFx(Istr-2:Iend+2,Jstr-2:Jend+2),  DVnew,
     &          VFe(Istr-2:Iend+2,Jstr-2:Jend+2)
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
      real sustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real svstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real sustrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real svstrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /smsdat_sustrg/sustrg /smsdat_svstrg/svstrg
      real    sustrp(2), svstrp(2), sms_time(2)
      real    sms_cycle, sms_scale
      integer*4 itsms, sms_ncycle, sms_rec, lsusgrd,
     &        lsvsgrd,sms_tid, susid, svsid
       common /smsdat1/ sustrp, svstrp, sms_time
       common /smsdat2/ sms_cycle, sms_scale
       common /smsdat3/ itsms, sms_ncycle, sms_rec, lsusgrd,
     &                  lsvsgrd,sms_tid, susid, svsid
	real srflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
	common /forces_srflx/srflx
	 real slp(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
	 common /forces_slp/slp
         real slpg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
         common /slpdat_slpg/slpg
         real slpp(2), slp_time(2)
         real slp_cycle,slp_scale
         integer*4 itslp, slp_ncycle, slp_rec,
     &   lslpgrd, slp_tid, slp_id, slpunused
         common /slpdat/
     &  slpp, slp_time, slp_cycle, slp_scale,
     &  itslp, slp_ncycle, slp_rec, slp_tid,
     &  slp_id, lslpgrd,  slpunused
      real bustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real bvstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real bvstrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat/
     &        bms_tintrp, bustrp,       bvstrp,    tbms,
     &        bmsclen,    bms_tstart,   bms_tend,  tsbms,   sclbms,
     &        itbms,      bmstid,busid, bvsid,     tbmsindx,
     &        bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real visc2_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
     &       /mixing_visc2_sponge_p/visc2_sponge_p
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
      if (iic.eq.ntstart) then
        kbak=kstp
        kold=kstp
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      elseif (iic.eq.ntstart+1) then
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      else
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kstp-2
        if (kold.lt.1) kold=4
        cff1= 1.781105D0
        cff2=-1.06221D0
        cff3= 0.281105D0
      endif
      do j=JstrV-2,jend+1
        do i=IstrU-2,iend+1
          Drhs(i,j)=h(i,j)+cff1*zeta(i,j,kstp)+cff2*zeta(i,j,kbak)
     &                                         +cff3*zeta(i,j,kold)
        enddo
      enddo
      do j=jstr-1,jend+1
        do i=IstrU-1,iend+1
          urhs(i,j)=cff1*ubar(i,j,kstp) +cff2*ubar(i,j,kbak)
     &                                         +cff3*ubar(i,j,kold)
          DUon(i,j)=0.5D0*(Drhs(i,j)+Drhs(i-1,j))*on_u(i,j)*urhs(i,j)
        enddo
      enddo
      do j=JstrV-1,jend+1
        do i=istr-1,iend+1
          vrhs(i,j)=cff1*vbar(i,j,kstp) +cff2*vbar(i,j,kbak)
     &                                         +cff3*vbar(i,j,kold)
          DVom(i,j)=0.5D0*(Drhs(i,j)+Drhs(i,j-1))*om_v(i,j)*vrhs(i,j)
        enddo
      enddo
      if (iic.eq.ntstart) then
        cff0=1.D0
        cff1=0.D0
        cff2=0.D0
        cff3=0.D0
      elseif (iic.eq.ntstart+1) then
        cff0= 1.0833333333333D0
        cff1=-0.1666666666666D0
        cff2= 0.0833333333333D0
        cff3=0.D0
      else
        cff0=0.614D0
        cff1=0.285D0
        cff2=0.088D0
        cff3=0.013D0
      endif
      do j=JstrV-1,jend
        do i=IstrU-1,iend
          zeta_new(i,j)=zeta(i,j,kstp) + dtfast*pm(i,j)*pn(i,j)
     &            *(DUon(i,j)-DUon(i+1,j)+DVom(i,j)-DVom(i,j+1))
          zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
          Dnew(i,j)=zeta_new(i,j)+h(i,j)
          UFx(i,j)=cff0*zeta_new(i,j) +cff1*zeta(i,j,kstp)
     &             +cff2*zeta(i,j,kbak) +cff3*zeta(i,j,kold)
          UFe(i,j)=UFx(i,j)
          VFe(i,j)=UFx(i,j)*UFx(i,j)
        enddo
      enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=zeta_new(i,j)
        enddo
      enddo
      call zetabc_tile (Istr,Iend,Jstr,Jend)
       cff=0.5D0*g
        fac=0.5D0*100.0D0/rho0
      do j=jstr,jend
        do i=istr,iend
          rubar(i,j)=cff*on_u(i,j)*( (h(i-1,j)+h(i,j))*(UFe(i-1,j)
     &                        -UFe(i,j)) +VFe(i-1,j)-VFe(i,j)
     &                                                              )
          rvbar(i,j)=cff*om_v(i,j)*( (h(i,j-1)+h(i,j))*(UFe(i,j-1)
     &                        -UFe(i,j)) +VFe(i,j-1)-VFe(i,j)
     &                                                              )
           call get_slp
        rubar(i,j)=rubar(i,j)+ fac*on_u(i,j)* (h(i-1,j)+h(i,j)+
     &             UFe(i-1,j)+UFe(i,j))*(slp(i-1,j)-slp(i,j))
        rvbar(i,j)=rvbar(i,j)+ fac*om_v(i,j)* (h(i,j-1)+h(i,j)+
     &             UFe(i,j-1)+UFe(i,j))*(slp(i,j-1)-slp(i,j))
        enddo
      enddo
      do j=jstr,jend
        do i=istr-1,iend
          UFx(i,j)=0.25D0*(DUon(i,j)+DUon(i+1,j))
     &                     *(urhs(i,j)+urhs(i+1,j))
          VFx(i+1,j)=0.25D0*(DUon(i+1,j)+DUon(i+1,j-1))
     &                       *(vrhs(i+1,j)+vrhs(i,j))
     &                                 *pmask(i+1,j)
        enddo
      enddo
      do j=jstr-1,jend
        do i=istr,iend
          VFe(i,j)=0.25D0*(DVom(i,j)+DVom(i,j+1))
     &                      *(vrhs(i,j)+vrhs(i,j+1))
          UFe(i,j+1)=0.25D0*(DVom(i,j+1)+DVom(i-1,j+1))
     &                       *(urhs(i,j+1)+urhs(i,j))
     &                                 *pmask(i,j+1)
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          rubar(i,j)=rubar(i,j)-UFx(i,j)+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i,j)
          rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j)+VFe(i,j-1)
        enddo
      enddo
      do j=JstrV-1,jend
        do i=IstrU-1,iend
          cff=Drhs(i,j)*(
     &                             fomn(i,j)
     &  +0.5D0*( dndx(i,j)*(vrhs(i,j)+vrhs(i,j+1))
     &        -dmde(i,j)*(urhs(i,j)+urhs(i+1,j)))
     &                                          )
          UFx(i,j)=cff*(vrhs(i,j)+vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i,j)+urhs(i+1,j))
        enddo
      enddo
      do j=jstr,jend
        do i=IstrU,iend
          rubar(i,j)=rubar(i,j)+0.25D0*(UFx(i,j)+UFx(i-1,j))
        enddo
      enddo
      do j=JstrV,jend
        do i=istr,iend
          rvbar(i,j)=rvbar(i,j)-0.25D0*(VFe(i,j)+VFe(i,j-1))
        enddo
      enddo
      do j=jstr-1,jend
        do i=istr-1,iend
          cff=2.D0*Drhs(i,j)*visc2_r(i,j)
          UFx(i,j)=cff*(ubar(i+1,j,kstp)-ubar(i,j,kstp))
     &                                 *pm(i,j)*on_r(i,j)
          VFe(i,j)=cff*(vbar(i,j+1,kstp)-vbar(i,j,kstp))
     &                                 *pn(i,j)*om_r(i,j)
          cff1=0.0625D0*visc2_p(i+1,j+1)*( Drhs(i,j)
     &       +Drhs(i+1,j)+Drhs(i,j+1)+Drhs(i+1,j+1) )*(
     &          (pn(i+1,j+1)+pn(i,j+1)+pn(i+1,j)+pn(i,j))
     &             *(ubar(i+1,j+1,kstp)-ubar(i+1,j,kstp))
     &         +(pm(i+1,j+1)+pm(i,j+1)+pm(i+1,j)+pm(i,j))
     &             *(vbar(i+1,j+1,kstp)-vbar(i,j+1,kstp))
     &                                                  )
     &                     *pmask(i+1,j+1)
          UFe(i+1,j+1)=cff1*om_p(i+1,j+1)
          VFx(i+1,j+1)=cff1*on_p(i+1,j+1)
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          rubar(i,j)=rubar(i,j)+UFx(i,j)-UFx(i-1,j)
     &                         +UFe(i,j+1)-UFe(i,j)
          rvbar(i,j)=rvbar(i,j)+VFx(i+1,j)-VFx(i,j)
     &                         +VFe(i,j)-VFe(i,j-1)
        enddo
      enddo
      if (rdrg2.gt.0.D0) then
        do j=jstr,jend
          do i=IstrU,iend
            cff=0.25D0*( vbar(i  ,j,kstp)+vbar(i  ,j+1,kstp)
     &                +vbar(i-1,j,kstp)+vbar(i-1,j+1,kstp))
            rubar(i,j)=rubar(i,j)-ubar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,kstp)*ubar(i,j,kstp)+cff*cff)
     &                                     )*om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,jend
          do i=istr,iend
            cff=0.25D0*( ubar(i,j  ,kstp)+ubar(i+1,j  ,kstp)
     &                +ubar(i,j-1,kstp)+ubar(i+1,j-1,kstp))
            rvbar(i,j)=rvbar(i,j)-vbar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,kstp)*vbar(i,j,kstp))
     &                                     )*om_v(i,j)*on_v(i,j)
          enddo
        enddo
      else if (rdrg.gt.0.0D0) then
        do j=jstr,jend
          do i=IstrU,iend
            rubar(i,j)=rubar(i,j) - rdrg*ubar(i,j,kstp)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,jend
          do i=istr,iend
            rvbar(i,j)=rvbar(i,j) - rdrg*vbar(i,j,kstp)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
      endif
      do j=JstrV-1,jend
        do i=IstrU-1,iend
          DUon(i,j)=zeta(i,j,kstp)+h(i,j)
        enddo
      enddo
      cff=0.5D0*dtfast
      cff2=2.D0*dtfast
      do j=jstr,jend
        do i=IstrU,iend
          DUnew=( (DUon(i,j)+DUon(i-1,j))*ubar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
     &                          *rubar(i,j)+cff2*sustr(i,j)
     &                                                    )
     &                                         *umask(i,j)
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
        enddo
      enddo
      do j=JstrV,jend
        do i=istr,iend
          DVnew=( (DUon(i,j)+DUon(i,j-1))*vbar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
     &                          *rvbar(i,j)+cff2*svstr(i,j)
     &                                                    )
     &                                         *vmask(i,j)
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
        enddo
      enddo
      call    u2dbc_tile (istr,iend,jstr,jend, UFx)
      call    v2dbc_tile (istr,iend,jstr,jend, UFx)
      call diag_tile (istr,iend,jstr,jend, UFx,UFe)
      call exchange_r2d_tile (istr,iend,jstr,jend,
     &                   zeta(-1,-1,knew))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                   ubar(-1,-1,knew))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                   vbar(-1,-1,knew))
      VMAXL=100.D0
      VMAX=0.D0
      do j=JstrV,Jend
        do i=Istr,Iend
          IF(ABS(vbar(i,j,knew)).GE.VMAX) THEN
            VMAX=ABS(vbar(i,j,knew))
            IMAX=i+iminmpi-1
            JMAX=j+jminmpi-1
          ENDIF
        enddo
      enddo
      IF(VMAX.GT.VMAXL) THEN
         write(stdout,'(A)') '                                         '
         write(stdout,'(A)') '                                         '
         write(stdout,'(A)') ' ======================================= '
         write(stdout,'(A)') ' =                                     = '
         write(stdout,'(A)') ' =   STEP2D:   ABNORMAL JOB END        = '
         write(stdout,'(A)') ' =                 BLOW UP             = '
         write(stdout,'(A)') ' =                                     = '
         write(stdout,'(A)') ' ======================================= '
         write(stdout,'(A)') '                                         '
         write(stdout,'(A,F8.2)') '  VMAX (M/S) =',VMAX
         write(stdout,'(A,2I6)')  '  IMAX JMAX  =',IMAX,JMAX
        WRITE(stdout,'(A,I6)')    '  IIC        =',iic
        may_day_flag=1
        stop
      ENDIF
      return
      end
