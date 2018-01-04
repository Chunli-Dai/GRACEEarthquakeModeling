        module global
        implicit none
!!      INTEGER,save    :: NPT
!!      real*8,save,allocatable ::LAT(:),LON(:),R(:),alph(:),bet(:)
!!      REAL*8,save, ALLOCATABLE:: sgg_LORF(:,:)
!       IN JPL GSM RL05 headline: 
!       EARTH 3.9860044150e+14 6.3781363000e+06
        real*8,parameter:: R0 = 6378136.3d0     !Earth mean radius [m]
        real*8,parameter:: GM = 398600.4415d9   !Grav. const * Earth mass [m3s-2]
        real*8,save :: PI
!!      PI= 4.d0*datan(1.d0)
        end module global
