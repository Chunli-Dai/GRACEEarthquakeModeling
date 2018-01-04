c layer one and two flipped, after the read statement!
c layer 1: water
c layer 2: ice

      parameter(ityp=360)
      parameter(nla=90,nlo=180)

      dimension fvel(ityp,8),fvels(ityp,8),frho(ityp,8)
      dimension fthi(ityp,7)
      dimension amapvp(8,nlo,nla),amaprho(8,nlo,nla),
     +          amapvs(8,nlo,nla),amapthi(7,nlo,nla),
     +          amapele(nlo,nla)
      character*2 ctype(ityp),line*506,dum*1,dum0*5
      character*2 types(nlo),nsta*4,ntyp*2,atype(nlo,nla)
      character*12 names(7)
	double precision latN,latS,lonW,lonE
        integer ilatN,ilatS,ilonW,ilonE
        double precision thi_ave(8),vp_ave(8),vs_ave(8),rho_ave(8)
        double precision elev_ave
        integer cnt,err
      data names/'water','ice','soft sed.','hard sed.',
     +         'upper crust','middle crust','lower crust'/
      open(2,file='CNtype2_key.txt')
      open(7,file='CNtype2.txt')
      open(8,file='CNelevatio2.txt')

      dx=360/nlo
c... read in key for crust types
c...............................
      read(2,890)dum
      print*,' ... reading key file ...'
      do 101 i=1,ityp
         read(2,899)ctype(i)
c        print 899,ctype(i)
         read(2,899)line
         read(line,*)(fvel(i,l),l=1,8)
         read(2,899)line
         read(line,*)(fvels(i,l),l=1,8)
         read(2,899)line
         read(line,*)(frho(i,l),l=1,8)
         read(2,899)line
         read(line,*)(fthi(i,l),l=1,7)
c flip layers
         aux=fvel(i,1)
         fvel(i,1)=fvel(i,2)
         fvel(i,2)=aux
         aux=fvels(i,1)
         fvels(i,1)=fvels(i,2)
         fvels(i,2)=aux
         aux=frho(i,1)
         frho(i,1)=frho(i,2)
         frho(i,2)=aux
         aux=fthi(i,1)
         fthi(i,1)=fthi(i,2)
         fthi(i,2)=aux
 101  continue

c... read CNtype file
c...............................
      read(7,*)flons
      print*,' ... reading model ...'
      read(8,899)line
      do 40 j=1,nla
         read(8,*)ilat,(amapele(i,j),i=1,nlo)
         read(7,901)ilat,types
         do 10 i=1,nlo
            do 20 l=1,ityp
            if(types(i).eq.ctype(l))then
              atype(i,j)=ctype(l)
              do 30 k=1,8
              amapvp(k,i,j)=fvel(l,k)
              amapvs(k,i,j)=fvels(l,k)
              amaprho(k,i,j)=frho(l,k)
 30           continue
              do 31 k=1,7
 31           amapthi(k,i,j)=fthi(l,k)
              goto 10
            endif
 20         continue
            print*,' crust type code not found: ',types(i)
            print*,' latitude: ',ilat,' long index: ',i
 10      continue
 40   continue

*-------------------
c     now look up coordinates
     
      open(66,file='outcr') 
      print*,' the output file is outcr'
      print*,' '
 60   continue
!     print*,' enter lat, lon  (* quits)'
!     read(*,'(a)')line
!     if(line(1:1).eq.' '.or.line(1:1).eq.'*')goto 99
!     read(line,*,err=99)flat,flon
!     cola=90.-flat
!     if(flon.gt.180.)flon=flon-360.
!     ilat=int(cola/dx)+1
!     ilon=int((flon+180.)/dx)+1

!     flatarry=(/34,36,38,40,42/)
!     flonarry=(/138,140,142,144,146/)      

!     write(*,*)'Enter the latN latS lonW lonE'
!     read(*,*)latN,latS,lonW,lonE
!     if(lonW.gt.180.)lonW=lonW-180.d0
!     if(lonE.gt.180.)lonE=lonE-180.d0
!     ilatN=int((90d0-latN)/dx)+1;ilatS=int((90d0-latS)/dx)+1+1;
!     ilonW=int((lonW+180.)/dx)+1;ilonE=int((lonE+180.)/dx)+1+1;

       open(1,file='elev_select.dat')
       cnt=0
       do while(.true.)
        read(1,*,iostat=err)
        if (err.lt.0)exit
        cnt=cnt+1
       enddo
       write(*,*)'Number of selected points:',cnt
       rewind(1)
!     do ilat=ilatN,ilatS
!       flat=90d0-dble(ilat-1)*dx
!     do ilon= ilonW,ilonE
!       flon=dble(ilon-1)*dx-180d0
      thi_ave=0d0;vp_ave=0d0;vs_ave=0d0;rho_ave=0d0
      elev_ave=0d0
      do j=1,cnt
        read(1,*)flon,flat
      cola=90.-flat
      if(flon.gt.180.)flon=flon-360.
      ilat=int(cola/dx)+1
      ilon=int((flon+180.)/dx)+1
      print 999,ilat,ilon,atype(ilon,ilat)
 999  format(2i5,1x,a2)
      vthi=0.
      vvp=0.
      vvs=0.
      vrho=0.
      do 50 i=2,7
         vthi=vthi+amapthi(i,ilon,ilat)
         vvp=vvp+amapthi(i,ilon,ilat)/amapvp(i,ilon,ilat)
         vvs=vvs+amapthi(i,ilon,ilat)/amapvs(i,ilon,ilat)
         vrho=vrho+amapthi(i,ilon,ilat)*amaprho(i,ilon,ilat)
 50   continue
      vvp=vthi/vvp
      vvs=vthi/vvs
      vrho=vrho/vthi
      write(66,793)'type, latitude, longitude, elevation: ',
     +     atype(ilon,ilat),flat,flon,amapele(ilon,ilat)
      write(66,794)'crustal thickness, ave. vp, vs, rho:  ',
     +     '  ',vthi,vvp,vvs,vrho                                
      write(66,793)'Mantle below Moho: ave. vp, vs, rho:  ',
     +     '  ',amapvp(8,ilon,ilat),amapvs(8,ilon,ilat),
     +     amaprho(8,ilon,ilat)
 793  format(a,2x,a2,2x,11x,4f11.4)
 794  format(a,2x,a2,2x,4f11.4)
      write(66,*)' ' 
      write(66,'(a)')' 7-layer crustal model (thickness, vp,vs,rho)'
      if(amapele(ilon,ilat).lt.0..and.amapthi(1,ilon,ilat).ne.0)   
     +  amapthi(1,ilon,ilat)=-amapele(ilon,ilat)/1000.
      do 70 i=1,7
         write(66,792)amapthi(i,ilon,ilat),amapvp(i,ilon,ilat),
     +        amapvs(i,ilon,ilat),amaprho(i,ilon,ilat),names(i)
 70   continue
 791  format(a2,2x,2f11.4,7f6.2)
 792  format(4f10.4,2x,a12)
!     goto 60
      elev_ave=elev_ave+amapele(ilon,ilat)
      do i=1,7
        thi_ave(i)=thi_ave(i)+amapthi(i,ilon,ilat)
        vp_ave(i)=vp_ave(i)+amapthi(i,ilon,ilat)/amapvp(i,ilon,ilat)
        vs_ave(i)=vs_ave(i)+amapthi(i,ilon,ilat)/amapvs(i,ilon,ilat)
        rho_ave(i)=rho_ave(i)+amapthi(i,ilon,ilat)*amaprho(i,ilon,ilat)
      enddo
!       write(*,*)'vs water=',vs_ave(1) !Infinity
      do i=8,8
        vp_ave(i)=vp_ave(i)+1d0/amapvp(i,ilon,ilat)
        vs_ave(i)=vs_ave(i)+1d0/amapvs(i,ilon,ilat)
        rho_ave(i)=rho_ave(i)+amaprho(i,ilon,ilat)
      enddo
!     write(1,'(3f11.4)')flon,flat,amapele(ilon,ilat)*1d-3  !meter to km
!     enddo !ilon
!     enddo !ilat
      enddo !j
      
      elev_ave=elev_ave/dble(cnt)
      do i=1,7
        thi_ave(i)=thi_ave(i)/dble(cnt)
        if(vp_ave(i).gt.1d-13)then !
          vp_ave(i)=thi_ave(i)/(vp_ave(i)/dble(cnt))
        else !in case all thickness for this layer is zero
          vp_ave(i)=amapvp(i,ilon,ilat) !Choose the last pixel as the average
        endif
        if(vs_ave(i).gt.1d-13)then
          vs_ave(i)=thi_ave(i)/(vs_ave(i)/dble(cnt))
        else
          vs_ave(i)=amapvs(i,ilon,ilat)
        endif
        if(thi_ave(i).gt.1d-13)then
          rho_ave(i)=rho_ave(i)/dble(cnt)/thi_ave(i)
        else
          rho_ave(i)=amaprho(i,ilon,ilat)
        endif
      enddo
!       write(*,*)'vs water=',vs_ave(1)
      do i=8,8
        vp_ave(i)=1d0/(vp_ave(i)/dble(cnt))
        vs_ave(i)=1d0/(vs_ave(i)/dble(cnt))
        rho_ave(i)=rho_ave(i)/dble(cnt)
      enddo
      close(1)
      !For both ocean and land area, thi_ave(1) consider the ocean at land as zero
      !Use the average elevation instead
      if(elev_ave.lt.0..and.thi_ave(1).ne.0) then 
      !Consider Bason area: elev_ave <0, thi_ave(1)==0
        thi_ave(1)=-elev_ave/1000.d0
      elseif(elev_ave.ge.0)then
      !Works for coastal area. Not for mountain lake
      !Mountain lake: elev_ave>0, but thi_ave(1) should not be zero
      !Ignore this mountain lake case
        thi_ave(1)=0d0
      endif

      open(1,file='crust2_select_ave_weight.dat')
!     write(1,*)' 7-layer crustal model (thickness, vp,vs,rho); Weighted
!    + Average over the 2011 Tohoku EQ area'
      write(1,'(a,f11.4)')' 7-layer crustal model (thickness,vp,vs,rho)
     +; Elevation (m): ',elev_ave
      do i=1,7
        write(1,792)thi_ave(i),vp_ave(i),vs_ave(i),rho_ave(i),names(i)
      enddo
      i=8
      write(1,'(10x,3f10.4,2x,a12)')vp_ave(i),vs_ave(i),rho_ave(i),
     +'Mantle below Moho'
 890  format(////a)
 899  format(a)
 901  format(i4,1x,180(2x,a2,1x))

 99   continue
      close(7)
      close(66)
	close(1)

      end
