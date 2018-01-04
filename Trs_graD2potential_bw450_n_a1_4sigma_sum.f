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
!! Modified from GNR_STATC_fast_JP_var_gra.f by dcl Sep 20, 2012
!! Purpose: truncate the gravity signal from the seismic model,
!!          ~/Osuwork/GRACE/Seismic2gravityByLei/WangLei20100515_dai/
!!          GravityChangeofFaultIntegrationFitGRACE/FaultIntegration/caltech_graD_spherical.dat
!!          to spherical harmonic degree 60.
!!        Transform it using f(θ,λ)=sumsum[Clmcos(mλ)+Slmsin(mλ)]Plm(cosθ) to
!!          match the Grid2SHCs.exe 
!! Modified from Seismic_gra_truncate_nmax60.f by dcl Sep 20, 2012
!! Purpose: Recover the gravity signal from the truncated spherical harmonic
!!          coefficients
!! Modified from Seismic_gra_truncate_nmax60_SHCs2Grid.f by dcl, Oct 2, 2012
!! Purpose: Transform the coefficients expanded from gravity Down component
!!          :caltech_graD_spherical_NMAX60.coe
!!          to Standard Potential coefficients:
!!          caltech_graD_spherical_NMAX60_potential.coe
!! Modified from Trs_graD2potential_bw450.f
!! Purpose: to transform n files together
!!
!!          


CC	use module global

	IMPLICIT NONE
	INTEGER,ALLOCATABLE:: NLON(:),IS(:),IE(:),IS2(:),IE2(:),IS3(:)
     +				,IE3(:),INUM(:)
	REAL*8, ALLOCATABLE:: PBAR(:),PBAR1(:),PBAR2(:)
	REAL*8, ALLOCATABLE:: IBDN(:),ATPY(:),EST(:),RES(:),
     +				DRC(:),QVEC(:),SVEC(:)
	INTEGER            :: NMAX,SITERATION,NITERATION,NCOEFS,NCOMP2
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
	REAL*8, ALLOCATABLE:: ATPA(:,:),J2N(:)
	INTEGER            :: IA,JA,lk
        double precision,allocatable ::cbar(:),sbar(:)
        double precision,allocatable ::cbar2(:),sbar2(:)

!Aug 9, 2012
!       parameter(NMAX=60)
!       parameter(NMAX=449)
        CHARACTER*120 :: IDIR,ODIR,filenames,ifilet
        CHARACTER*120 :: ifilet2
        character(len=80),allocatable :: ifile(:),ofile(:)
        character(len=80),allocatable :: IFILEsigma(:),ofilesigma(:)
        character(len=80),allocatable :: ofilesum(:)
        integer inm,num_model
!       NAMELIST /parm/ IDIR,ODIR,filenames
        INTEGER,ALLOCATABLE::IS1(:),IE1(:)
        INTEGER NCOEFS1
        REAL*8  graxyz(3),covgra(3)
        REAL*8, ALLOCATABLE:: DISP(:)
        double precision,allocatable ::cstd(:),sstd(:)
        double precision,allocatable ::cstd2(:),sstd2(:)
        REAL*8, ALLOCATABLE:: EST2(:),DISP2(:)
        integer npt
        REAL*8  D2R,TPI
        real*8 timex,time0
!Aug 29, 2012
        character gridfile*80
!       NAMELIST /parm/ IDIR,ODIR,filenames,gridfile
!       NAMELIST /parm/ NMAX,IDIR,ODIR,filenames
!Sep 20, 2012
        real*8,allocatable :: gra(:),APT(:,:)
        real*8 graest(1)
!Oct 2, 2012
        real*8, parameter :: Re=6371d3  !The Earth Radius 
! in ~/Osuwork/GRACE/Seismic2gravityByLei/WangLei20100515_dai/
!    GravityChangeofFaultIntegrationFitGRACE/FaultIntegration/GravityInt.f
        integer N
        real*8 :: cnm, snm
        real*8 :: cnm2, snm2
        real*8 oceanthi  !unit: meter        
        real*8 a1
        NAMELIST /parm/ NMAX,IDIR,ODIR,filenames,oceanthi
!

        PI= 4.d0*datan(1.d0); 
        D2R=PI/180D0;TPI=2d0*PI
!
! EITHER STANDARD OR NAMELIST INPUT
!
        READ(unit=*, NML=parm, IOSTAT=ios) 
        a1=Re-oceanthi;
        write(*,*)'Radius of Ocean bottom:',a1,'meter'
        
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
        ALLOCATE(cbar2(NCOEFS1),sbar2(NCOEFS1))
	ALLOCATE(EST2(NCOEFS),DISP2(NCOEFS))
        ALLOCATE(cstd(NCOEFS1),sstd(NCOEFS1))
        ALLOCATE(cstd2(NCOEFS1),sstd2(NCOEFS1))
!       ALLOCATE(ATPA(NCOEFS,NCOEFS),ATPY(NCOEFS))
	
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

! Open coefficient file
        
        OPEN (UNIT=1,FILE=filenames,STATUS='OLD')
        read(1,*)num_model
        ALLOCATE(IFILE(num_model),OFILE(num_model))
        ALLOCATE(IFILEsigma(num_model),OFILEsigma(num_model))
        ALLOCATE(OFILEsum(num_model))
        do inm=1,num_model
          read(1,*)IFILE(inm),OFILE(inm)
          read(1,*)IFILEsigma(inm),OFILEsigma(inm)
          read(1,*)OFILEsum(inm)
        
          ifilet=trim(IDIR)//trim(IFILE(inm))
          ifilet2=trim(IDIR)//trim(IFILEsigma(inm))
          write(*,*)'inm=',inm
          write(*,*)'coe file=',ifilet

!         ifilet='caltech_graD_spherical_NMAX60.coe'
!         ifilet='caltech_graD_spherical_bw450.coe'
!         write(*,*)'coe file=',ifilet
        
!         call read_GOCO(ifilet,nmax,NCOEFS1,cbar,sbar)
          call read_SHstd(ifilet,nmax,NCOEFS1,cbar,sbar,cstd,sstd)
          call read_SHstd(ifilet2,nmax,NCOEFS1,cbar2,sbar2,cstd2,sstd2)
       
!!        Set zero- and first-degree coefficients to zeros
!         cbar(1)=0d0;cbar(2)=0d0;
!         cbar(NMAX+2)=0d0;sbar(NMAX+2)=0d0;
!         cstd(1)=0d0;cstd(2)=0d0;
!         cstd(NMAX+2)=0d0;sstd(NMAX+2)=0d0;
 
!       Assign cbar sbar to EST
           EST=0D0
           DISP=0D0
           EST2=0D0
           DISP2=0D0
           M=0
           L=NMAX-M+1
           DO K=1,L
               EST(IS(M+1)+K-1) =cbar(IS1(m+1)+K-1)
               DISP(IS(M+1)+K-1) =cstd(IS1(m+1)+K-1)

               EST2(IS(M+1)+K-1) =cbar2(IS1(m+1)+K-1)
               DISP2(IS(M+1)+K-1) =cstd2(IS1(m+1)+K-1)
           enddo
           DO M=1,NMAX
             L=NMAX-M+1
             DO K=1,L
               EST(IS(M+1)+K-1) =cbar(IS1(m+1)+K-1)
               EST(IS(M+1)+L+K-1)=sbar(IS1(m+1)+K-1)
               DISP(IS(M+1)+K-1) =cstd(IS1(m+1)+K-1)
               DISP(IS(M+1)+L+K-1)=sstd(IS1(m+1)+K-1)

               EST2(IS(M+1)+K-1) =cbar2(IS1(m+1)+K-1)
               EST2(IS(M+1)+L+K-1)=sbar2(IS1(m+1)+K-1)
               DISP2(IS(M+1)+K-1) =cstd2(IS1(m+1)+K-1)
               DISP2(IS(M+1)+L+K-1)=sstd2(IS1(m+1)+K-1)
             enddo
           enddo
!          DISP=0D0
 
!        write(*,*)'Check whether assign the COE correctly:'
! 	M=0;
! 	DO K=IS3(M+1),IE3(M+1)
! 	WRITE(*,'(2I5,4(1x,E23.15))')K-IS3(M+1)+M,M,EST(K),0D0
!     +          ,DISP(K),0d0
! 	ENDDO
! 	DO M=1,NMAX
! 	DO K=IS3(M+1),IS3(M+1)+INUM(M+1)-1
! 	WRITE(*,'(2I5,4(1x,E23.15))')K-IS3(M+1)+M,M,EST(K)
!     +          ,EST(K+INUM(M+1))
!     +          ,DISP(K),DISP(K+INUM(M+1))
! 	ENDDO
! 	ENDDO
!!
!        write(*,*)'Check whether assign the COEsigma correctly:'
! 	M=0;
! 	DO K=IS3(M+1),IE3(M+1)
! 	WRITE(*,'(2I5,4(1x,E23.15))')K-IS3(M+1)+M,M,EST2(K),0D0
!     +          ,DISP2(K),0d0
! 	ENDDO
! 	DO M=1,NMAX
! 	DO K=IS3(M+1),IS3(M+1)+INUM(M+1)-1
! 	WRITE(*,'(2I5,4(1x,E23.15))')K-IS3(M+1)+M,M,EST2(K)
!     +          ,EST2(K+INUM(M+1))
!     +          ,DISP2(K),DISP2(K+INUM(M+1))
! 	ENDDO
! 	ENDDO
!!
         OPEN(2,file=trim(ODIR)//trim(OFILE(inm)),status='unknown')
         OPEN(3,file=trim(ODIR)//trim(OFILEsigma(inm)),status='unknown')
         open(8,file=trim(ODIR)//trim(OFILEsum(inm)),status='unknown')
	WRITE(*,*)'OUTPUT factor'
!       open(2,file='caltech_graD_spherical_NMAX60_potential.coe')
!       open(2,file='caltech_graD_spherical_bw450_potential.coe')
        M=0;
        DO K=IS3(M+1),IE3(M+1)
          N=K-IS3(M+1)+M
!       WRITE(2,'(2I5,2E23.15)') K-IS3(M+1)+M,M,EST(K),0D0
        !First from microGal to m/s^2
          cnm=EST(K)*1d-8*(a1/R0)**(N+2)/dble(N+1)/(GM/R0/R0)
          cnm2=EST2(K)*(a1/R0)**N*(a1/Re)**2
!	write(*,*)n,m,(a1/R0)**N*(a1/Re)**2
        WRITE(2,'(2I5,2E23.15)') N,M,cnm,0D0
        WRITE(3,'(2I5,2E23.15)') N,M,cnm2,0D0
        WRITE(8,'(2I5,2E23.15)') N,M,cnm+cnm2,0D0
        ENDDO
        DO M=1,NMAX
        DO K=IS3(M+1),IS3(M+1)+INUM(M+1)-1
          N=K-IS3(M+1)+M
          cnm=EST(K)*1d-8*(a1/R0)**(N+2)/dble(N+1)/(GM/R0/R0)
          snm=EST(K+INUM(M+1))*1d-8*(a1/R0)**(N+2)/dble(N+1)/(GM/R0/R0)
          cnm2=EST2(K)*(a1/R0)**N*(a1/Re)**2
          snm2=EST2(K+INUM(M+1))*(a1/R0)**N*(a1/Re)**2
!	write(*,*)n,m,(a1/R0)**N*(a1/Re)**2
        WRITE(2,'(2I5,2E23.15)') N,M,cnm,snm
        WRITE(3,'(2I5,2E23.15)') N,M,cnm2,snm2
        WRITE(8,'(2I5,2E23.15)') N,M,cnm+cnm2,snm+snm2
        ENDDO
        ENDDO
        close(2)
        close(3)
        close(8)
        enddo ! do inm=1,num_model
        close(1)
        DEALLOCATE(IFILE,OFILE)

        call timexe(timex); timex=timex-time0
        print*, "execution time [sec]:",real(timex)

	STOP
	END PROGRAM 

  	Include "lgdr2.f90"
        include "brmul.f"
        include 'brinv.f'
        include 'read_SHstd.f'


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

