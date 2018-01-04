      subroutine subAPgrad(APgrad,NMAX,npt,LATin,LONin,R)
        use global
        implicit none
        integer M,L,i,npt,PT,NMAX,NCOEFS
        integer lk,K,ii,jj
        INTEGER,ALLOCATABLE:: IS(:),IE(:)
        REAL*8, ALLOCATABLE:: PBAR(:),PBAR1(:),PBAR2(:)
        REAL*8, ALLOCATABLE:: AP(:,:) !,APgrad(:,:,:)
!       REAL*8, ALLOCATABLE:: LATin(:),LONin(:)
        REAL*8, ALLOCATABLE:: LAT(:),LON(:)
        REAL*8  D2R
        REAL*8 latin(npt),lonin(npt),R(npt)
        REAL*8 APgrad(6,(NMAX+1)**2,npt)

        PI= 4.d0*datan(1.d0);
        D2R=PI/180D0;

        NCOEFS = (NMAX+1)**2
        ALLOCATE(AP(9,NCOEFS))  !,APgrad(6,NCOEFS,npt))
        allocate(lat(npt),lon(npt))  !,latin(npt),lonin(npt))
        ALLOCATE(IS(NMAX+1),IE(NMAX+1))

        do i=1,npt
           lat(i)=90d0-latin(i);
           lat(i)=lat(i)*D2R; lon(i)=lonin(i)*D2R
        enddo

!       R=R0

        IS(1)=1
        IE(1)=NMAX+1
        DO I=2,NMAX+1
           IS(I)=IE(I-1)+1
           IE(I)=IS(I)+2*(NMAX-I+2)-1
        ENDDO

        do i=1,npt
          PT=I
          AP=0d0
          include "AP.f"

!         do jj=1,NCOEFS
            !!AP:xx,xy,xz,yx,yy,yz,zx,zy,zz
            !!APgrad: XX,yy,zz,xy,xz,yz
            APgrad(1,:,i)=AP(1,:)
            APgrad(2,:,i)=AP(5,:)
            APgrad(3,:,i)=AP(9,:)
            APgrad(4,:,i)=AP(2,:)
            APgrad(5,:,i)=AP(3,:)
            APgrad(6,:,i)=AP(6,:)
!         enddo

        enddo !i= npt
        return
      end

