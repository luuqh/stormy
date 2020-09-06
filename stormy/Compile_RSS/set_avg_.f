      subroutine set_avg (tile)
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
      integer*4 tile, trd
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
      call set_avg_tile (Istr,Iend,Jstr,Jend)
      return
      end
      subroutine set_avg_tile (Istr,Iend,Jstr,Jend)
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
      integer*4 Istr,Iend,Jstr,Jend, i,j, ilc
      real cff
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
      real zeta_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real ubar_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real vbar_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_zeta/zeta_avg
     &       /avg_ubar/ubar_avg
     &       /avg_vbar/vbar_avg
      real bostr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_bostr/bostr_avg
      real wstr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_wstr/wstr_avg
      real sustr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_sustr/sustr_avg
      real svstr_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_svstr/svstr_avg
      real visc2_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
     &       /mixing_visc2_sponge_p/visc2_sponge_p
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
      ilc=1+iic-ntstart
      if (ilc.gt.ntsavg) then
        if (mod(ilc-1,navg).eq.1) then
          if (wrtavg(indxZ)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                zeta_avg(i,j)=zeta(i,j,kstp)
              enddo
            enddo
          endif
          if (wrtavg(indxUb)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                ubar_avg(i,j)=ubar(i,j,kstp)
              enddo
            enddo
          endif
        elseif (mod(ilc-1,navg).gt.1) then
          if (wrtavg(indxZ)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                zeta_avg(i,j)=zeta_avg(i,j)+zeta(i,j,kstp)
              enddo
            enddo
          endif
          if (wrtavg(indxUb)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                ubar_avg(i,j)=ubar_avg(i,j)+ubar(i,j,kstp)
              enddo
            enddo
          endif
          if (wrtavg(indxVb)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                vbar_avg(i,j)=vbar_avg(i,j)+vbar(i,j,kstp)
              enddo
            enddo
          endif
          if (wrtavg(indxBostr)) then
            do j=Jstr,Jend
              do i=Istr,Iend
                bostr_avg(i,j)=bostr_avg(i,j)+
     &                         0.5D0*sqrt((bustr(i,j)+bustr(i+1,j))**2
     &                                 +(bvstr(i,j)+bvstr(i,j+1))**2)
     &                                                          *rho0
              enddo
            enddo
          endif
          if (wrtavg(indxWstr)) then
            do j=Jstr,Jend
              do i=Istr,Iend
                wstr_avg(i,j)=wstr_avg(i,j)+
     &                         0.5D0*sqrt((sustr(i,j)+sustr(i+1,j))**2
     &                                 +(svstr(i,j)+svstr(i,j+1))**2)
     &                                                          *rho0
              enddo
            enddo
          endif
          if (wrtavg(indxUWstr)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                sustr_avg(i,j)=sustr_avg(i,j)+
     &                         sustr(i,j)*rho0
              enddo
            enddo
          endif
          if (wrtavg(indxVWstr)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                svstr_avg(i,j)=svstr_avg(i,j)+
     &                         svstr(i,j)*rho0
              enddo
            enddo
          endif
        elseif (mod(ilc-1,navg).eq.0) then
          cff=1.D0/float(navg)
          if (Istr+Jstr.eq.2) time_avg=time_avg+float(navg)*dt
          if (wrtavg(indxZ)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                zeta_avg(i,j)=cff*( zeta_avg(i,j)
     &                   +zeta(i,j,kstp))
              enddo
            enddo
          endif
          if (wrtavg(indxUb)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                ubar_avg(i,j)=cff*( ubar_avg(i,j)
     &                   +ubar(i,j,kstp))
              enddo
            enddo
          endif
          if (wrtavg(indxVb)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                vbar_avg(i,j)=cff*( vbar_avg(i,j)
     &                   +vbar(i,j,kstp))
              enddo
            enddo
          endif
          if (wrtavg(indxBostr)) then
            do j=Jstr,Jend
              do i=Istr,Iend
                bostr_avg(i,j)=cff*( bostr_avg(i,j)+
     &                         0.5D0*sqrt((bustr(i,j)+bustr(i+1,j))**2
     &                                 +(bvstr(i,j)+bvstr(i,j+1))**2)
     &                                                          *rho0)
              enddo
            enddo
          endif
          if (wrtavg(indxWstr)) then
            do j=Jstr,Jend
              do i=Istr,Iend
                wstr_avg(i,j)=cff*( wstr_avg(i,j)+
     &                         0.5D0*sqrt((sustr(i,j)+sustr(i+1,j))**2
     &                                 +(svstr(i,j)+svstr(i,j+1))**2)
     &                                                          *rho0)
              enddo
            enddo
          endif
          if (wrtavg(indxUWstr)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                sustr_avg(i,j)=cff*(sustr_avg(i,j)+
     &                         sustr(i,j)*rho0)
              enddo
            enddo
          endif
          if (wrtavg(indxVWstr)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                svstr_avg(i,j)=cff*(svstr_avg(i,j)+
     &                         svstr(i,j)*rho0)
              enddo
            enddo
          endif
        endif
      endif
      return
      end
