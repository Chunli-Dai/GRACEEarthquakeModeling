Program Ggrid
implicit none
integer i,j
double precision lat,lon
integer nmax,bw
double precision dB,dL
double precision latN,latS,lonW,lonE
!Modified Feb 5, 2013
!input bw and gridfile
character gridfile*80

!bw=240
write(*,*)'Input bw, gridfile'
read(*,*)bw,gridfile
write(*,*)'Input latN latS lonW lonE:'
read(*,*)latN,latS,lonW,lonE

dB=180d0/2d0/dble(bw)
dL=2d0*dB

!open(1,file='grid')
open(1,file=gridfile)
!!geodetic latitude, longitude
!do i=1,180
!  do j=1,360
!    lat=-89.5d0+dble(i-1)
!    lon=0.5d0+dble(j-1)
!    write(1,fmt='(3F8.2)')lat,lon,0d0
!  enddo
!enddo

!do i=1,2*bw
!  do j=1,2*bw
!    lat=90.d0-dble(i-.5d0)*dB
!    lon=dble(j-5d0/8d0)*dL
!    write(1,fmt='(3F10.4)')lat,lon,0d0
!  enddo
!enddo

do i=1,2*bw
  do j=1,2*bw
    lat=90.d0-dble(i-.5d0)*dB
!   lon=dble(j-.5d0)*dL
    lon=dble(j-1d0)*dL
    if(lat.ge.latS.and.lat.le.latN.and.lon.ge.lonW.and.lon.le.lonE)then
    write(1,fmt='(3F10.4)')lat,lon,0d0
    endif
  enddo
enddo


close(1)

stop
end
