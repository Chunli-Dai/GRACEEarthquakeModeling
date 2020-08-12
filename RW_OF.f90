        Program RW_OF
        !Purpose: read OF_bw900.txt
        !         write it with format
        integer i,j,signland,err
        double precision lat, lon
        character ifile*80,ofile*80
        double precision elev
        double precision latN,latS,lonW,lonE
        double precision dx

        write(*,*)'Input latN latS lonW lonE:'
        read(*,*)latN,latS,lonW,lonE

        if(lonW.lt.0)lonW=lonW+360d0
        if(lonE.lt.0)lonE=lonE+360d0

        !Expand the region by 5 degrees.
        dx=5d0;
        latN=latN+dx; latS=latS-dx; lonW=lonW-dx;lonE=lonE+dx;

!       ifile='etopo1.xyz';ofile='etopo1_JPfmt.xyz'
        ifile='ETOPO1_Bed_g_int.xyz';ofile='etopo1_JPfmt.xyz'

        open(1,file=ifile)
        open(2,file=ofile)
        do while(.true.)
!       read(1,fmt=*,iostat=err)lat,lon,signland
        read(1,fmt=*,iostat=err)lon,lat,elev
!164.783333333333331 15.0166666666666657 -5226
        if(err.lt.0)exit
        if(lon.lt.0)lon=lon+360d0
        if(lat.ge.latS.and.lat.le.latN.and.lon.ge.lonW.and.lon.le.lonE)then
!         write(2,fmt='(f6.2,2x,f6.2,2x,i1)')lat,lon,signland
          write(2,fmt='(f20.15,2x,f20.15,2x,f7.0)')lat,lon,elev
        endif
        enddo

        close(1)
        close(2)
        stop
        end
