	PROGRAM GRACE_L2_PCG
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
!! Aug 10,2012
!! Modified from /0/home/dai.56/share/CHL/progdcl/GNR_STATC_fast_JP_var.f
!! Purpose: subtract a reference field and get the corresponding std

CC	use module global

	IMPLICIT NONE
	INTEGER,ALLOCATABLE:: NLON(:),IS(:),IE(:),IS2(:),IE2(:),IS3(:)
     +				,IE3(:),INUM(:)
	INTEGER            :: NMAX,SITERATION,NITERATION,NCOEFS,NCOMP2
	INTEGER		   :: I,J,K,L,M,ITR,CNT,PT,LATN,NLAT
	REAL*8             :: DLAT,PLAT,INC,P1ST,P2ND,P3RD,PA,PB,PC
	REAL*8             :: ALPHA,BETA,DELTA_OLD,DELTA_NEW,QFACTOR
	CHARACTER(LEN=1)   :: CH(2)
! ADDED BY YIQUN CHEN
	CHARACTER*120      :: inputfile
	INTEGER            :: ios,sign
	INTEGER            :: IA,JA,lk
        double precision,allocatable ::cbar(:),sbar(:)

!Aug 9, 2012
!       parameter(NMAX=60)
!       parameter(NMAX=120)
        CHARACTER*120 :: IDIR,ODIR,filenames,ifilet
        character(len=80) ifile,ofile
        integer inm,num_model
!       NAMELIST /parm/ IDIR,ODIR,filenames,NMAX
        integer ncut
        character C20file*160
	integer ostdflag,C20flag,C20rightflag
        character label*8
        NAMELIST /parm/ IDIR,ODIR,filenames,NMAX,ncut,C20file,ostdflag
     +          ,C20flag
        INTEGER,ALLOCATABLE::IS1(:),IE1(:)
        INTEGER NCOEFS1
        integer npt
        REAL*8  D2R,TPI
        real*8 timex,time0
        double precision,allocatable :: cbar2(:),sbar2(:)
        double precision,allocatable :: cstd(:),sstd(:)
        double precision,allocatable :: cstd2(:),sstd2(:)
        double precision dc,ds,dcstd,dsstd
        integer n
        
        character string*90
        integer err,nC20
        double precision,allocatable ::year(:),C20a(:),dC20(:),stdC20(:)
        double precision mjdt,yeart,C20t,dC20t,stdC20t
        integer yeari,doyi
        logical alive

!
! EITHER STANDARD OR NAMELIST INPUT
!
        READ(unit=*, NML=parm, IOSTAT=ios) 
        
        write(*,*)'idir=',IDIR,ODIR,filenames,NMAX,ncut,C20file,ostdflag
     +          ,C20flag
! INITIALIZATION
!
        NCOEFS=(nmax+2)*(nmax+1)/2
        ALLOCATE(IS1(NMAX+1),IE1(NMAX+1),cbar(NCOEFS),sbar(NCOEFS))
	ALLOCATE(IS(NMAX+1),IE(NMAX+1))
        ALLOCATE(cbar2(NCOEFS),sbar2(NCOEFS))
        ALLOCATE(cstd(NCOEFS),sstd(NCOEFS))
        ALLOCATE(cstd2(NCOEFS),sstd2(NCOEFS))
	
!
! DEFINE INDICES
!
!
        IS(1)=1
        IE(1)=NMAX+1
        DO I=2,NMAX+1
          IS(I)=IE(I-1)+1
          IE(I)=IS(I)+(NMAX-I+2)-1
        ENDDO

! Read C20 file
      if(C20flag.eq.1)then
        inquire(file=C20file,exist=alive)
        if(.not.alive)then
          write(*,*)C20file," doesn't exist."
          stop
        endif

        open(1,file=C20file,STATUS='OLD')
        read(1,fmt='(a90)')string
        i=0
        do while(.true.)
          read(1,fmt='(a80)',iostat=err)string
          if (err.lt.0)exit
          i=i+1   
        enddo
        nC20=i
        write(*,*)'Number of C20:',nC20
        allocate(year(nC20),C20a(nC20),dC20(nC20),stdC20(nC20))
        rewind(1)
        read(1,fmt='(a90)')string
        do i=1,nC20
          read(1,fmt='(a80)',iostat=err)string
          if (err.lt.0)exit
            read(string,fmt=*,iostat=err)mjdt,yeart,C20t,dC20t,stdC20t
            year(i)=yeart;C20a(i)=C20t;dC20(i)=dC20t*1.0E-10;
            stdC20(i)=stdC20t*1.0E-10
        enddo
        close(1)
      endif

!
! Open coefficient file
        OPEN (UNIT=1,FILE=filenames,STATUS='OLD')
        read(1,*)label,IFILE
        ifilet=trim(IDIR)//trim(IFILE)
        call read_SHstd(label,ifilet,nmax,NCOEFS,cbar,sbar,cstd,sstd)
          do m=0,nmax
            do k=IS(m+1),IE(m+1)
               n=m+k-IS(m+1)
        write(*,*)n,m,cbar(k),sbar(k),cstd(k),sstd(k)
            enddo
         enddo

        read(1,*)num_model
        do inm=1,num_model
          read(1,*)label,IFILE,OFILE
          !2003_016_0.txt  2003_016_2_M.txt
          if(C20flag.eq.1)then
          read(IFILE,'(i4,1x,i3)')yeari,doyi
          yeart=dble(yeari)+dble(doyi-15)/365.25d0;
!         if(C20flag.eq.1)then
            if(abs(year(inm)-yeart).gt.10./365.25d0)then
              write(*,*)'Wrong C20 for:',IFILE
!             C20rightflag=0 !Bad idea to use the CSR solutions since the
!             background C20 is different from the mean value in C20 file
!             Better to use the value month before.
!             Avoid this case
              C20rightflag=1
            else
              C20rightflag=1
            dC20t=dC20(inm)
            stdC20t=stdC20(inm)
            endif
          endif

          ifilet=trim(IDIR)//trim(IFILE)
          write(*,*)'inm=',inm
          write(*,*)'coe file=',ifilet
        
       call read_SHstd(label,ifilet,nmax,NCOEFS,cbar2,sbar2,cstd2,sstd2)

           OPEN(2,file=trim(ODIR)//trim(OFILE),status='unknown')
 
           do m=0,nmax
             do k=IS(m+1),IE(m+1)
                n=m+k-IS(m+1)
!        write(*,*)n,m,cbar2(k),sbar2(k),cstd2(k),sstd2(k)
                dc=cbar2(k)-cbar(k);ds=sbar2(k)-sbar(k);
                dcstd=dsqrt(cstd2(k)**2+cstd(k)**2)
                dsstd=dsqrt(sstd2(k)**2+sstd(k)**2)
                if(N.lt.ncut)then
                  dc=0d0; ds=0d0
                  dcstd=0d0;dsstd=0d0
                endif
                if(C20flag.eq.1.and.C20rightflag.eq.1)then
                if(N.eq.2.and.M.eq.0)then
                  dc=dC20t;dcstd=stdC20t
                  ds=0d0;dsstd=0d0
                endif
                elseif(C20flag.eq.2)then
                if(N.eq.2.and.M.eq.0)then
                  dc=0d0;dcstd=0d0
                  ds=0d0;dsstd=0d0
                endif
                endif
                if(ostdflag.eq.1)then
                  write(2,fmt='(2I5,4(1x,E23.15))')n,m,dc,ds,dcstd,dsstd
		elseif(ostdflag.eq.0)then
		  write(2,fmt='(2I5,2(1x,E23.15))')n,m,dc,ds
	        endif
             enddo
           enddo
 
         close(2)

!!Another output order
!        open(3,file=trim(ODIR)//trim(OFILE),status='unknown')
!        do n=0,NMAX
!          do m=0,n
!            k=IS(m+1)+n-m
!        write(*,*)n,m,cbar2(k),sbar2(k),cstd2(k),sstd2(k)
!            dc=cbar2(k)-cbar(k);ds=sbar2(k)-sbar(k);
!            dcstd=dsqrt(cstd2(k)**2+cstd(k)**2)
!            dsstd=dsqrt(sstd2(k)**2+sstd(k)**2)
!             if(N.lt.ncut)then
!               dc=0d0; ds=0d0
!               dcstd=0d0;dsstd=0d0
!             endif
!            write(3,fmt='(2I5,4(1x,E23.15))')n,m,dc,ds,dcstd,dsstd
!          enddo
!        enddo
!        close(3)


        enddo ! do inm=1,num_model
        close(1)
	
	STOP
	END PROGRAM 

        include 'read_SHstd_wwoH.f'

