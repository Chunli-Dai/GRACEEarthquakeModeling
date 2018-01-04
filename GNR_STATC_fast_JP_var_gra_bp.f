        include "CNST_GOCE.f"
	PROGRAM GRACE_L2_PCG
	use global
CC Modified from SOFT.GRACE_L2.PCG.f90 by dcl-2010-1-29
CC Purpose: to estimate parameters using Conjugate Gradient Method
CC The annotations "CC" are made by dcl

CC Modified from GOCE_INC.f by dcl-2010-2-13
CC Modified from GOCE_TEST.f by dcl-2010-2-18
CC Modification:	call gnr_statc_sub.f90
CC			instead of read sggobs_test
!! Modified from test/GOCE_TEST3_woab6.f by dcl Aug 9, 2012
!! Purpose: compared with GNR_STATC_fast_JP.f90, compute
!!          the Gradients tensor from GRACE L2 RL05 data for JP region;
!!          compute the variance for the Gradients tensor;
!! Modified from GNR_STATC_fast_JP_var.f by dcl Aug 29, 2012
!! Purpose: compute the gravity vector in NED frame from GRACE L2 RL05 data
!!          compute the variance for the three gravity components
!!Modified from GNR_STATC_fast_JP_var_gra_woset0_bw450.f , dcl, Feb 5, 2013
!!Purpose: add nmax to parm
!!Modified from GNR_STATC_fast_JP_var_gra_woset0_nmax.f , dcl, June 6, 2014
!!Purpose: add all cases to one, compatible to older versions
!!Modification: 1\ add ncut, rc, grid fmt type to parm;
!!              2\ Use read_SHstd_wwoH.f instead for different SH format
!!Version 2, Nov, 2014
!!Modification: 1\ store A matrix first for faster computation

CC	use module global

	IMPLICIT NONE
	INTEGER,ALLOCATABLE:: NLON(:),IS(:),IE(:),IS2(:),IE2(:),IS3(:)
     +				,IE3(:),INUM(:)
	REAL*8, ALLOCATABLE:: PBAR(:),PBAR1(:),PBAR2(:)
CC	REAL*8, ALLOCATABLE:: PNM(:),PNM1(:),PNM2(:),PNMC(:),QNM(:),QNM1(:),QNM2(:),QNMC(:)
CC	REAL*8, ALLOCATABLE:: LAT1(:),LON1(:),R1(:),LAT2(:),LON2(:),R2(:),MLAT1(:),MLAT2(:)
	REAL*8, ALLOCATABLE:: IBDN(:),ATPY(:),EST(:),RES(:),
     +				DRC(:),QVEC(:),SVEC(:)
	INTEGER            :: NMAX,SITERATION,NITERATION,NCOEFS,NCOMP2
CC	INTEGER            :: I,J,K,L,M,ITR,CNT,NPT,PT,LATN,NLAT
	INTEGER		   :: I,J,K,L,M,ITR,CNT,PT,LATN,NLAT
	REAL*8             :: DLAT,PLAT,INC,P1ST,P2ND,P3RD,PA,PB,PC
	REAL*8             :: ALPHA,BETA,DELTA_OLD,DELTA_NEW,QFACTOR
	CHARACTER(LEN=1)   :: CH(2)
! ADDED BY YIQUN CHEN
	CHARACTER*120      :: inputfile
	INTEGER            :: ios,sign
  	REAL*8, ALLOCATABLE:: LAT(:),LON(:),R(:),alph(:),bet(:)
	REAL*8, ALLOCATABLE:: AP(:,:),ASRP(:,:),ASRPT(:,:),ASR(:,:)
	REAL*8, ALLOCATABLE:: ATPAT(:,:)
CC	REAL*8, ALLOCATABLE:: data_LORF(:,:)
	REAL*8, ALLOCATABLE:: ATPA(:,:),J2N(:)
	INTEGER            :: IA,JA,lk
        double precision,allocatable ::cbar(:),sbar(:)

!Aug 9, 2012
!       parameter(NMAX=60)
!       parameter(NMAX=449)
        CHARACTER*120 :: IDIR,ODIR,filenames,ifilet
        character(len=80),allocatable :: ifile(:),ofile(:)
        integer inm,num_model
!       NAMELIST /parm/ IDIR,ODIR,filenames
        INTEGER,ALLOCATABLE::IS1(:),IE1(:)
        INTEGER NCOEFS1
        REAL*8  graxyz(3),covgra(3)
        REAL*8, ALLOCATABLE:: DISP(:)
        double precision,allocatable ::cstd(:),sstd(:)
        integer npt
        REAL*8  D2R,TPI
        real*8 timex,time0
!Aug 29, 2012
        character gridfile*80,label*8
        integer ncut,N,grdflag,err
        real*8 rc
        REAL*8, ALLOCATABLE:: APgra(:,:,:)
        
!       NAMELIST /parm/ IDIR,ODIR,filenames,gridfile
       NAMELIST /parm/ NMAX,IDIR,ODIR,filenames,gridfile,ncut,rc,grdflag

        PI= 4.d0*datan(1.d0); 
        D2R=PI/180D0;TPI=2d0*PI
!
! EITHER STANDARD OR NAMELIST INPUT
!
        READ(unit=*, NML=parm, IOSTAT=ios) 
        
        write(*,*)'Inputed information, NMAX=',NMAX
        write(*,*)'IDIR:',IDIR
        write(*,*)'ODIR:',ODIR
        write(*,*)'ncut:',ncut
        write(*,*)'rc:',rc, ' meter'
        write(*,*)'grdflag:',grdflag

        if (grdflag.eq.0)then
        open(unit=1,file=gridfile,STATUS="OLD")
        cnt=0
        do while(.true.)
                read(1,*,iostat=err); if (err.lt.0)exit
                cnt=cnt+1
        enddo
        npt=cnt
        allocate(lat(npt),lon(npt),r(npt))
        rewind(1)
        do i=1,npt
                read(1,*) lat(i),lon(i)
        !       read(1,*) lat(i),lon(i); r(i)=R0+250d3          !dcl, use the upper
        !       line, 2009-11-9
                lat(i)=90d0-lat(i); lat(i)=lat(i)*D2R; lon(i)=lon(i)*D2R
        enddo
        elseif (grdflag.eq.1)then !format provide by JSG
          open(unit=1,file=gridfile,STATUS="OLD")
        read(1,*) !Header
        cnt=0
        do while(.true.)
                read(1,*,iostat=err); if (err.lt.0)exit
                cnt=cnt+1
        enddo
        npt=cnt
        allocate(lat(npt),lon(npt),r(npt))
        rewind(1)
        read(1,*) !Header
        do i=1,npt
                read(1,*) lon(i),lat(i)
                lat(i)=90d0-lat(i); lat(i)=lat(i)*D2R; lon(i)=lon(i)*D2R
        enddo
        endif
!        r=R0
        close(1)
        print*, npt," POINTS READ."

        if(abs(rc).gt.100)then
          r=rc
        else
          r=R0
        endif
        
!
! INITIALIZATION
!
!       NCOEFS1 IS1 IE1 is from NCOEFS IS IE in GNR_STATC_fast_JP.f90
        NCOEFS1=(nmax+2)*(nmax+1)/2
        ALLOCATE(IS1(NMAX+1),IE1(NMAX+1),cbar(NCOEFS1),sbar(NCOEFS1))
	NCOEFS = (NMAX+1)**2
	NCOMP2 = nint(4D0/3D0*NMAX**3 + 3D0*NMAX**2+8D0/3D0*NMAX+1D0)
	ALLOCATE(IS(NMAX+1),IE(NMAX+1),IS2(NMAX+1),IE2(NMAX+1),
     +			IS3(NMAX+1),IE3(NMAX+1),INUM(NMAX+1))
	ALLOCATE(EST(NCOEFS))
        ALLOCATE(DISP(NCOEFS))
        ALLOCATE(cstd(NCOEFS1),sstd(NCOEFS1))
	
!
! DEFINE INDICES
!
!       FOR NORMAL MATRIX
!
        IS(1)=1
        IE(1)=NMAX+1
        DO I=2,NMAX+1
                IS(I)=IE(I-1)+1
                IE(I)=IS(I)+2*(NMAX-I+2)-1
        ENDDO
!
!       FOR INVERSE OF B.D. NORMAL MATRIX, IBDN 
!
        IS2(1)=1
        IE2(1)=(NMAX+1)**2
        DO I=2,NMAX+1
                IS2(I)=IE2(I-1)+1
                IE2(I)=IS2(I)+(2*(NMAX-I+2))**2-1
        ENDDO
!
!       FOR WEIGHTED OBSERVATION VECTOR, ATPY   
!
        IS3(1)=1
        IE3(1)=NMAX+1
        DO I=2,NMAX+1
                IS3(I)=IE3(I-1)+1
                IE3(I)=IS3(I)+2*(NMAX-I+2)-1
        ENDDO
!

        IS1(1)=1
        IE1(1)=NMAX+1
        DO I=2,NMAX+1
          IS1(I)=IE1(I-1)+1
          IE1(I)=IS1(I)+(NMAX-I+2)-1
        ENDDO

!       FOR PARAMETER VECTOR, EST       
!
        INUM(1)=NMAX+1
        DO I=2,NMAX+1
                INUM(I)=NMAX-I+2
        ENDDO

!
          call timexe(time0)

!       'Prepare for gra and grad...'
        do i=1,npt
           lat(i)=lat(i)/D2R; lon(i)=lon(i)/D2R
           lat(i)=90d0-lat(i);
        enddo
        ALLOCATE(APgra(3,NCOEFS,npt))
        call subAPgra(APgra,NMAX,npt,LAT,LON,R)
        do i=1,npt
           lat(i)=90d0-lat(i); lat(i)=lat(i)*D2R; lon(i)=lon(i)*D2R
        enddo

! Open coefficient file
        OPEN (UNIT=1,FILE=filenames,STATUS='OLD')
        read(1,*)num_model
        ALLOCATE(IFILE(num_model),OFILE(num_model))
        do inm=1,num_model
          read(1,*)label,IFILE(inm),OFILE(inm)
        
          ifilet=trim(IDIR)//trim(IFILE(inm))
          write(*,*)'inm=',inm
          write(*,*)'coe file=',ifilet
        
!         call read_GOCO(ifilet,nmax,NCOEFS1,cbar,sbar)
          call read_SHstd(label,ifilet,nmax,NCOEFS1,cbar,sbar,cstd,sstd)
       
!!        Set zero- and first-degree coefficients to zeros
!         cbar(1)=0d0;cbar(2)=0d0;
!         cbar(NMAX+2)=0d0;sbar(NMAX+2)=0d0;
!         cstd(1)=0d0;cstd(2)=0d0;
!         cstd(NMAX+2)=0d0;sstd(NMAX+2)=0d0;
 
!       Assign cbar sbar to EST
          EST=0D0
          DISP=0D0
          M=0
          L=NMAX-M+1
          DO K=1,L
              N=K+M-1
              EST(IS(M+1)+K-1) =cbar(IS1(m+1)+K-1)
              DISP(IS(M+1)+K-1) =cstd(IS1(m+1)+K-1)
              if(N.lt.ncut)then
                EST(IS(M+1)+K-1) =0d0;DISP(IS(M+1)+K-1) =0d0
              endif
          enddo
          DO M=1,NMAX
            L=NMAX-M+1
            DO K=1,L
              N=K+M-1
              EST(IS(M+1)+K-1) =cbar(IS1(m+1)+K-1)
              EST(IS(M+1)+L+K-1)=sbar(IS1(m+1)+K-1)
              DISP(IS(M+1)+K-1) =cstd(IS1(m+1)+K-1)
              DISP(IS(M+1)+L+K-1)=sstd(IS1(m+1)+K-1)
              if(N.lt.ncut)then
                EST(IS(M+1)+K-1) =0d0;DISP(IS(M+1)+K-1) =0d0
                EST(IS(M+1)+L+K-1)=0d0;DISP(IS(M+1)+L+K-1)=0d0
              endif
            enddo
          enddo

        write(*,*)'Check whether assign the COE correctly:'
	M=0;
	DO K=IS3(M+1),IE3(M+1)
	WRITE(*,'(2I5,4(1x,E23.15))')K-IS3(M+1)+M,M,EST(K),0D0
     +          ,DISP(K),0d0
	ENDDO
	DO M=1,NMAX
	DO K=IS3(M+1),IS3(M+1)+INUM(M+1)-1
	WRITE(*,'(2I5,4(1x,E23.15))')K-IS3(M+1)+M,M,EST(K)
     +          ,EST(K+INUM(M+1))
     +          ,DISP(K),DISP(K+INUM(M+1))
	ENDDO
	ENDDO

          OPEN(2,file=trim(ODIR)//trim(OFILE(inm)),status='unknown')

          do i=1,npt
	PT=I
        ALLOCATE(AP(3,NCOEFS))
!	include "AP_gra.f"
        AP=APgra(:,:,I)
        call BRMUL(AP,EST,3,NCOEFS,1,graxyz)
!       call BRMUL(ASRP,DISP,6,NCOEFS,NCOEFS,AMP)
!       call BRMUL(AMP,ASRPT,6,NCOEFS,6,covgra)
        covgra=0d0
        do j=1,3
        do k=1,NCOEFS
          covgra(j)=covgra(j)+AP(j,k)**2*DISP(k)**2
        enddo
        enddo

!	North East Down
!       write(2,'(2F10.4,12(E17.8E2))')90d0-lat(i)/D2R,lon(i)/D2R
!    +          ,graxyz(1:6)*1d12,dsqrt(covgra(1:6))*1d12
!1d8 m/s^2->microgals
        write(2,'(2F10.4,6(E17.8E2))')90d0-lat(i)/D2R,lon(i)/D2R
     +          ,graxyz(1:3)*1d8,dsqrt(covgra(1:3))*1d8 

CC	DEALLOCATE(AP,ASRP,ASRPT,ASR)
	DEALLOCATE(AP)	
	
	ENDDO ! do i=1,npt
        close(2)
        enddo ! do inm=1,num_model
        close(1)
        DEALLOCATE(IFILE,OFILE)
	
	DEALLOCATE(LAT,LON,R)

        call timexe(timex); timex=timex-time0
        print*, "execution time [sec]:",real(timex)

	STOP
	END PROGRAM 

  	Include "lgdr2.f90"
        include "brmul.f"
!       include 'read_SHstd.f'
        include 'read_SHstd_wwoH.f'
        include 'subAPgra.f'


CC	SUBROUTINE INT2CHAR(DIGIT,INT,CH)
CC	IMPLICIT NONE
CC	INTEGER :: INT,DUM,DIGIT,I
CC	CHARACTER(LEN=1) :: CH(DIGIT)
CC	DUM=INT
CC	DO I=1,DIGIT
CC		CH(DIGIT+1-I)=CHAR(MODULO(INT,10**I)/10**(I-1)+48)
CC		INT=INT-MODULO(INT,10**I)
CC	ENDDO
CC	INT=DUM
CC	RETURN
CC	END SUBROUTINE INT2CHAR

        SUBROUTINE TIMEXE(TIME)
        REAL*8 TIME
        REAL*4 T(2)
        DATA TOT/0.D0/
        TIME=DTIME(T)+TOT
        TOT=TIME
        RETURN
        END

