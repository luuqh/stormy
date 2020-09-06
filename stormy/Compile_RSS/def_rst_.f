      subroutine def_rst (ncid, total_rec, ierr)
      implicit none
      logical create_new_file
      integer*4 ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  u2dgrd(3), v2dgrd(3),  auxil(2),  checkdims
     &      , b3dgrd(4)
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
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 512)
      parameter (nf_max_attrs = 4096)
      parameter (nf_max_vars = 4096)
      parameter (nf_max_name = 128)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 32)
      parameter(maxncdim = 100)
      parameter(maxncatt = 2000)
      parameter(maxncvar = 2000)
      parameter(maxncnam = 128)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      ierr=0
      lstr=lenstr(rstname)
      if (nrpfrst.gt.0) then
        lvar=total_rec - (1+mod(total_rec-1, nrpfrst))
        call insert_time_index (rstname, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefhis
      if (ncid.ne.-1) create_new_file=.false.
      if (mynode.gt.0) create_new_file=.false.
  10  if (create_new_file) then
        ierr=nf_create(rstname(1:lstr), nf_clobber, ncid)
        if (ierr.ne.nf_noerr) then
          write(stdout,'(/3(1x,A)/)') 'ERROR in DEF_RST: Cannot',
     &             'create restart NetCDF file:', rstname(1:lstr)
          goto 99
        endif
        if (nrpfrst.eq.0) total_rec=0
        call put_global_atts (ncid, ierr)
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,  r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        b3dgrd(1)=r2dgrd(1)
        b3dgrd(2)=r2dgrd(2)
        b3dgrd(4)=timedim
        if (total_rec.eq.1) call def_grid (ncid, r2dgrd)
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 rstTstep)
        ierr=nf_put_att_text (ncid, rstTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                              NF_DOUBLE, 1, timedim, rstTime)
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, rstTime, 'long_name', lvar,
     &                                  vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, rstTime, 'units',     lvar,
     &                                  vname(3,indxTime)(1:lvar))
        lvar=lenstr (vname(4,indxTime))
        ierr=nf_put_att_text(ncid, rstTime, 'field',     lvar,
     &                                  vname(4,indxTime)(1:lvar))
        lvar=lenstr(vname(1,indxZ))
        ierr=nf_def_var (ncid, vname(1,indxZ)(1:lvar),
     &                                  NF_DOUBLE, 3, r2dgrd,rstZ)
        lvar=lenstr(vname(2,indxZ))
        ierr=nf_put_att_text (ncid, rstZ, 'long_name', lvar,
     &                                     vname(2,indxZ)(1:lvar))
        lvar=lenstr(vname(3,indxZ))
        ierr=nf_put_att_text (ncid, rstZ, 'units', lvar,
     &                                     vname(3,indxZ)(1:lvar))
        lvar=lenstr(vname(4,indxZ))
        ierr=nf_put_att_text (ncid, rstZ, 'field', lvar,
     &                                     vname(4,indxZ)(1:lvar))
        lvar=lenstr(vname(1,indxUb))
        ierr=nf_def_var (ncid, vname(1,indxUb)(1:lvar),
     &                               NF_DOUBLE, 3, u2dgrd, rstUb)
        lvar=lenstr(vname(2,indxUb))
        ierr=nf_put_att_text (ncid, rstUb, 'long_name', lvar,
     &                                   vname(2,indxUb)(1:lvar))
        lvar=lenstr(vname(3,indxUb))
        ierr=nf_put_att_text (ncid, rstUb, 'units', lvar,
     &                                   vname(3,indxUb)(1:lvar))
        lvar=lenstr(vname(4,indxUb))
        ierr=nf_put_att_text (ncid, rstUb, 'field', lvar,
     &                                   vname(4,indxUb)(1:lvar))
        lvar=lenstr(vname(1,indxVb))
        ierr=nf_def_var (ncid, vname(1,indxVb)(1:lvar),
     &                               NF_DOUBLE, 3, v2dgrd, rstVb)
        lvar=lenstr(vname(2,indxVb))
        ierr=nf_put_att_text (ncid, rstVb, 'long_name', lvar,
     &                                   vname(2,indxVb)(1:lvar))
        lvar=lenstr(vname(3,indxVb))
        ierr=nf_put_att_text (ncid, rstVb, 'units',     lvar,
     &                                   vname(3,indxVb)(1:lvar))
        lvar=lenstr(vname(4,indxVb))
        ierr=nf_put_att_text (ncid, rstVb, 'field',     lvar,
     &                                  vname(4,indxVb)(1:lvar))
        ierr=nf_enddef(ncid)
        write(*,'(6x,4A,1x,A,i4)') 'DEF_RST - Created new ',
     &              'netCDF file ''', rstname(1:lstr), '''.'
     &               ,' mynode =', mynode
      elseif (ncid.eq.-1) then
        ierr=nf_open (rstname(1:lstr), nf_write, ncid)
        if (ierr. eq. nf_noerr) then
          ierr=checkdims (ncid, rstname, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfrst.eq.0) then
              ierr=rec+1 - nrecrst
            else
              ierr=rec+1 - (1+mod(nrecrst-1, abs(nrpfrst)))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5/8x,A,I5,1x,A/)'
     &           ) 'DEF_RST WARNING: Actual number of records', rec,
     &             'in netCDF file',  '''',  rstname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpfrst.eq.0) then
              total_rec=rec+1
              if (mynode.gt.0) total_rec=total_rec-1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr. ne. nf_noerr) then
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
            write(stdout,'(/1x,4A, 1x,A,I4/)')     'DEF_RST ERROR: ',
     &     'Cannot open restart netCDF file ''',rstname(1:lstr),'''.'
     &               ,' mynode =', mynode
            goto 99
          endif
        endif
        ierr=nf_inq_varid (ncid, 'time_step', rstTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', rstname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), rstTime)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), rstname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxZ))
        ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), rstZ)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxZ)(1:lvar),  rstname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxUb))
        ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), rstUb)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxUb)(1:lvar), rstname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxVb))
        ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), rstVb)
        if (ierr. ne. nf_noerr) then
          write(stdout,1) vname(1,indxVb)(1:lvar), rstname(1:lstr)
          goto 99
        endif
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_RST -- Opened ',
     &             'existing restart file,  record =', rec
     &              ,' mynode =', mynode
      else
        ierr=nf_open (rstname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          write(stdout,'(/1x,4A, 1x,A,I4/)') 'DEF_RST ERROR: Cannot',
     &          'open restart netCDF file ''', rstname(1:lstr), '''.'
     &           ,' mynode =', mynode
          goto 99
        endif
      endif
   1  format(/1x,'DEF_RST ERROR: Cannot find variable ''',
     &               A, ''' in netCDF file ''', A, '''.'/)
      if (total_rec.eq.1) call wrt_grid (ncid, rstname, lstr)
  99  return
      end
