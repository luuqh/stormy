      subroutine MessPass2D_tile (Istr,Iend,Jstr,Jend, A)
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
      include 'mpif.h'
      real A(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer*4 Istr,Iend,Jstr,Jend, i,j, isize,jsize,   iter,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr, mdii,mdjj
      real buf_snd4(4),     buf_snd2(4),
     &     buf_rev4(4),     buf_rev2(4),
     &     buf_snd1(4),     buf_snd3(4),
     &     buf_rev1(4),     buf_rev3(4)
      integer*4 sub_X,size_X, sub_E,size_E
      parameter (sub_X=Lm,  size_X=7+2*sub_X,
     &           sub_E=Mm,  size_E=7+2*sub_E)
      real ibuf_sndN(0:size_X), ibuf_revN(0:size_X),
     &     ibuf_sndS(0:size_X), ibuf_revS(0:size_X),
     &     jbuf_sndW(0:size_E), jbuf_sndE(0:size_E),
     &     jbuf_revW(0:size_E), jbuf_revE(0:size_E)
      integer*4 imin,imax,ishft, jmin,jmax,jshft
      if (ii.eq.0 .and. Istr.eq.1) then
        imin=Istr-1
      else
        imin=Istr
      endif
      if (ii.eq.NP_XI-1 .and. Iend.eq.Lmmpi) then
        imax=Iend+1
      else
        imax=Iend
      endif
      ishft=imax-imin+1
      if (jj.eq.0 .and. Jstr.eq.1) then
        jmin=Jstr-1
      else
        jmin=Jstr
      endif
      if (jj.eq.NP_ETA-1 .and. Jend.eq.Mmmpi) then
        jmax=Jend+1
      else
        jmax=Jend
      endif
      jshft=jmax-jmin+1
      isize=2*ishft
      jsize=2*jshft
      do iter=0,1
       mdii=mod(ii+iter,2)
       mdjj=mod(jj+iter,2)
        if (mdii.eq.0) then
          if (WEST_INTER) then
            do j=jmin,jmax
              jbuf_sndW(j-jmin      )=A(1,j)
              jbuf_sndW(j-jmin+jshft)=A(2,j)
            enddo
            call MPI_Irecv (jbuf_revW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW, jsize, MPI_DOUBLE_PRECISION,
     &                       p_W, 1, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (EAST_INTER) then
            do j=jmin,jmax
              jbuf_sndE(j-jmin      )=A(Lmmpi-1,j)
              jbuf_sndE(j-jmin+jshft)=A(Lmmpi  ,j)
            enddo
            call MPI_Irecv (jbuf_revE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, 2, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdjj.eq.0) then
          if (SOUTH_INTER) then
            do i=imin,imax
              ibuf_sndS(i-imin      )=A(i,1)
              ibuf_sndS(i-imin+ishft)=A(i,2)
            enddo
            call MPI_Irecv (ibuf_revS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER) then
            do i=imin,imax
              ibuf_sndN(i-imin      )=A(i,Mmmpi-1)
              ibuf_sndN(i-imin+ishft)=A(i,Mmmpi  )
            enddo
            call MPI_Irecv (ibuf_revN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdii.eq.0) then
          if (SOUTH_INTER .and. WEST_INTER) then
            buf_snd1(1)=A(1,1)
            buf_snd1(2)=A(2,1)
            buf_snd1(3)=A(1,2)
            buf_snd1(4)=A(2,2)
            call MPI_Irecv (buf_rev1,4, MPI_DOUBLE_PRECISION,  p_SW,
     &                               6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1,4, MPI_DOUBLE_PRECISION,  p_SW,
     &                               5, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (NORTH_INTER .and. EAST_INTER) then
            buf_snd2(1)=A(Lmmpi-1,Mmmpi-1)
            buf_snd2(2)=A(Lmmpi  ,Mmmpi-1)
            buf_snd2(3)=A(Lmmpi-1,Mmmpi  )
            buf_snd2(4)=A(Lmmpi  ,Mmmpi  )
            call MPI_Irecv (buf_rev2,4, MPI_DOUBLE_PRECISION,  p_NE,
     &                               5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2,4, MPI_DOUBLE_PRECISION,  p_NE,
     &                               6, MPI_COMM_WORLD,        ierr)
          endif
        endif
        if (mdii.eq.1) then
          if (SOUTH_INTER .and. EAST_INTER) then
            buf_snd3(1)=A(Lmmpi-1,1)
            buf_snd3(2)=A(Lmmpi  ,1)
            buf_snd3(3)=A(Lmmpi-1,2)
            buf_snd3(4)=A(Lmmpi  ,2)
            call MPI_Irecv (buf_rev3,4, MPI_DOUBLE_PRECISION,  p_SE,
     &                               8, MPI_COMM_WORLD, req(7), ierr)
            call MPI_Send  (buf_snd3,4, MPI_DOUBLE_PRECISION,  p_SE,
     &                               7, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER .and. WEST_INTER) then
            buf_snd4(1)=A(1,Mmmpi-1)
            buf_snd4(2)=A(2,Mmmpi-1)
            buf_snd4(3)=A(1,Mmmpi  )
            buf_snd4(4)=A(2,Mmmpi  )
            call MPI_Irecv (buf_rev4, 4, MPI_DOUBLE_PRECISION, p_NW,
     &                               7, MPI_COMM_WORLD, req(8), ierr)
            call MPI_Send  (buf_snd4, 4, MPI_DOUBLE_PRECISION, p_NW,
     &                               8, MPI_COMM_WORLD,         ierr)
          endif
        endif
      enddo
      if (WEST_INTER) then
        call MPI_Wait (req(1),status(1,1),ierr)
        do j=jmin,jmax
          A(-1,j)=jbuf_revW(j-jmin)
          A( 0,j)=jbuf_revW(j-jmin+jshft)
        enddo
      endif
      if (EAST_INTER) then
        call MPI_Wait (req(2),status(1,2),ierr)
        do j=jmin,jmax
          A(Lmmpi+1,j)=jbuf_revE(j-jmin)
          A(Lmmpi+2,j)=jbuf_revE(j-jmin+jshft)
        enddo
      endif
      if (SOUTH_INTER) then
        call MPI_Wait (req(3),status(1,3),ierr)
        do i=imin,imax
          A(i,-1)=ibuf_revS(i-imin )
          A(i, 0)=ibuf_revS(i-imin+ishft)
        enddo
      endif
      if (NORTH_INTER) then
        call MPI_Wait (req(4),status(1,4),ierr)
        do i=imin,imax
          A(i,Mmmpi+1)=ibuf_revN(i-imin)
          A(i,Mmmpi+2)=ibuf_revN(i-imin+ishft)
        enddo
      endif
      if (SOUTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(5),status(1,5),ierr)
        A(-1,-1)=buf_rev1(1)
        A( 0,-1)=buf_rev1(2)
        A(-1, 0)=buf_rev1(3)
        A( 0, 0)=buf_rev1(4)
      endif
      if (NORTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(6),status(1,6),ierr)
        A(Lmmpi+1,Mmmpi+1)=buf_rev2(1)
        A(Lmmpi+2,Mmmpi+1)=buf_rev2(2)
        A(Lmmpi+1,Mmmpi+2)=buf_rev2(3)
        A(Lmmpi+2,Mmmpi+2)=buf_rev2(4)
      endif
      if (SOUTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(7),status(1,7),ierr)
        A(Lmmpi+1,-1)=buf_rev3(1)
        A(Lmmpi+2,-1)=buf_rev3(2)
        A(Lmmpi+1, 0)=buf_rev3(3)
        A(Lmmpi+2, 0)=buf_rev3(4)
      endif
      if (NORTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(8),status(1,8),ierr)
        A(-1,Mmmpi+1)=buf_rev4(1)
        A( 0,Mmmpi+1)=buf_rev4(2)
        A(-1,Mmmpi+2)=buf_rev4(3)
        A( 0,Mmmpi+2)=buf_rev4(4)
      endif
      return
      end
