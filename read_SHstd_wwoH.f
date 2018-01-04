	subroutine read_SHstd(label,ifile,nmax,NCOEFS,cbar,sbar,cstd,
     +                          sstd)
!!      Copy from Osuwork/classes/ES851 Sem Hydro & Oceano/
!!              Project/code/read_SH.f
!!	Purpose: read GOCO01S_gfc
!!               Read SH in this format for any order
!!                n m cnm snm
!!      Aug 10, 2012
!!      Modified from /0/home/dai.56/share/CHL/progdcl/read_SH2.f
!!      Purpose read Level2 GRACE data or  GOCO02S, as well as their sigmas
!!      Aug 13, 2012
!!      Modified from
!/0/home/dai.56/Osuwork/GRACE/Software/Subtract_Reference_SHCs/codedcl/read_SHstd.f
!!      Modification: no string before n m

	implicit none
	INTEGER :: nmax,i,j,l,m,NCOEFS,nmaxg,err,mm,k
	integer,allocatable :: IS(:),IE(:)
	double precision :: tmp,ct,st
!	double precision :: cbar((nmax+2)*(nmax+1)/2),sbar((nmax+2)*(nmax+1)/2)
!	double precision ,allocatable :: cbar(:),sbar(:)
	double precision cbar(NCOEFS),sbar(NCOEFS)
	character string*130,ifile*120
        logical alive
!!      Aug 10, 2012
        double precision cstd(NCOEFS),sstd(NCOEFS),cstdt,sstdt
        character ctype*8,label*8
        integer nlen
	

!	NCOEFS=(nmax+2)*(nmax+1)/2
	ALLOCATE(IS(NMAX+1),IE(NMAX+1))
!	ALLOCATE(cbar(NCOEFS),sbar(NCOEFS))

        IS(1)=1
        IE(1)=NMAX+1
        DO I=2,NMAX+1
                IS(I)=IE(I-1)+1
                IE(I)=IS(I)+(NMAX-I+2)-1
        ENDDO

	
!	open(1,file='2011_105_11m10_Deco.txt')
!       write(*,*)'ifile=',ifile
        inquire(file=ifile,exist=alive)
        if(.not.alive)then
          write(*,*)ifile," doesn't exist."
          stop
        endif

        open(101,file=ifile)
!	do while(.true.)
!	  read(1,fmt="(a80)",iostat=err)string
!	  if (err.lt.0)exit
!	  if (string(1:11).eq.'end_of_head')exit
!	enddo
	
	cbar=0d0;sbar=0d0
        cstd=0d0;sstd=0d0       !added Mar 9, 2013, 
                                !to avoid the wrong large std for
                                !2004_007_2_gra_std.txt

        nlen=len(trim(label))
        do while(.true.)
          read(101,fmt='(a130)',iostat=err)string
          if (err.lt.0)exit
          if (nlen.eq.0)then
            read(string,fmt=*,iostat=err)l,mm,ct,st,cstdt,sstdt
          elseif(string(1:nlen).eq.trim(label))then
            read(string,fmt=*,iostat=err)ctype,l,mm,ct,st,cstdt,sstdt
          else
            cycle
          endif
            if (l.ge.0.and.l.le.nmax)then
              K=IS(mm+1)+l-mm
              cbar(K)=ct
              sbar(K)=st         
              cstd(K)=cstdt
              sstd(K)=sstdt
            endif
        enddo
        close(101)

	return
	end
	
