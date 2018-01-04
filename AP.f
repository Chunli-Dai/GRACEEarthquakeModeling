        M=0
        L=NMAX+1
        ALLOCATE(PBAR(L),PBAR1(L),PBAR2(L))
        CALL LGDR2(DCOS(LAT(PT)),NMAX,M,PBAR,PBAR1,PBAR2)
        DO K=1,L
CC        PNMC(IS(1)+K-1) = (R0/R(PT))**(M+K) * (PBAR(K))
          lk=k+M-1
          AP(1,IS(1)+K-1) =1.D0/R(PT)**2*GM/R0*(R0/R(PT))**(lk+1)*
CC     +          (PBAR2(K)*sin(LAT(PT))**2
CC     +          -PBAR1(K)*cos(LAT(PT))-dble(lk+1)*PBAR(K))
CC use the next line, dcl-2010-2-13
     +		(PBAR2(K)-dble(lk+1)*PBAR(K))
          AP(2,IS(1)+K-1) =0.d0
          AP(3,IS(1)+K-1) =1.d0/R(PT)**2*GM/R0*(R0/R(PT))**(lk+1)
CC     +          *PBAR1(K)*sin(LAT(PT))*dble(-2-lk)
CC use the next line, dcl-2010-2-13
     +		*(-PBAR1(K))*dble(-2-lk)
          AP(4,IS(1)+K-1) =AP(2,IS(1)+K-1)
          AP(5,IS(1)+K-1) =-1.d0/R(PT)**2*GM/R0*(R0/R(PT))**(lk+1)*
CC     +          (dble(lk+1)*PBAR(K)+cos(LAT(PT))*PBAR1(K))
CC use the next line, dcl-2010-2-13
     +    (dble(lk+1)*PBAR(K)-Dcos(LAT(PT))*PBAR1(K)/Dsin(LAT(PT)))
          AP(6,IS(1)+K-1) =0.d0
          AP(7,IS(1)+K-1) =AP(3,IS(1)+K-1)
          AP(8,IS(1)+K-1) =AP(6,IS(1)+K-1)
          AP(9,IS(1)+K-1) =GM/R0**3*(R0/R(PT))**(lk+3)
     +           *dble((lk+1)*(lk+2))*PBAR(K)
CC	  AP(1,IS(1)+K-1) = -AP(5,IS(1)+K-1)-AP(9,IS(1)+K-1)
        ENDDO
        DEALLOCATE(PBAR,PBAR1,PBAR2)

        DO M=1,NMAX
        L=NMAX-M+1
        ALLOCATE(PBAR(L),PBAR1(L),PBAR2(L))
        CALL LGDR2(DCOS(LAT(PT)),NMAX,M,PBAR,PBAR1,PBAR2)
        DO K=1,L
CC      PNMC(IS(M+1)+K-1  ) = (R0/R(PT))**(M+K) * (PBAR(K)) * COS(M*LON(PT))
CC      PNMC(IS(M+1)+L+K-1) = (R0/R(PT))**(M+K) * (PBAR(K)) * SIN(M*LON(PT))
          lk=k+M-1      ! The number of l
          AP(1,IS(M+1)+K-1) =1.D0/R(PT)**2*GM/R0*(R0/R(PT))**(lk+1)
CC     +          *(PBAR2(K)*sin(LAT(PT))**2-PBAR1(K)*cos(LAT(PT))-
CC  use the next line,dcl-2010-2-13, found PBAR2 is the derivative relative to sita.
     +		*(PBAR2(K)-
     +          dble(lk+1)*PBAR(K))*DCOS(dble(M)*LON(PT))
          AP(1,IS(M+1)+L+K-1)=1.D0/R(PT)**2*GM/R0*(R0/R(PT))**(lk+1)
CC     +          *(PBAR2(K)*sin(LAT(PT))**2-PBAR1(K)*cos(LAT(PT))-
CC  use the next line,dcl-2010-2-13,
     +		*(PBAR2(K)-
     +          dble(lk+1)*PBAR(K))*DSIN(DBLE(M)*LON(PT))
          AP(2,IS(M+1)+K-1) = 1.d0/(R(PT)**2*Dsin(LAT(PT)))*GM/R0*
     +          (R0/R(PT))**(lk+1)*dble(m)*(PBAR(K)/Dtan(LAT(PT))+
CC     +          PBAR1(K)*sin(LAT(PT)))* SIN(M*LON(PT))
CC  use the next line,dcl-2010-2-13
     +		(-1.d0)*PBAR1(K))* DSIN(DBLE(M)*LON(PT))
          AP(2,IS(M+1)+L+K-1)=-1.d0/(R(PT)**2*Dsin(LAT(PT)))*GM/R0*
     +          (R0/R(PT))**(lk+1)*dble(m)*(PBAR(K)/Dtan(LAT(PT))+
CC     +          PBAR1(K)*sin(LAT(PT)))*COS(M*LON(PT))    
CC  use the next line,dcl-2010-2-13
     +		  (-1.d0)*PBAR1(K))*DCOS(DBLE(M)*LON(PT))
          AP(3,IS(M+1)+K-1) = 1.d0/R(PT)**2*GM/R0*
CC     +          (R0/R(PT))**(lk+1)*PBAR1(K)*sin(LAT(PT))
CC  use the next line,dcl-2010-2-13
     +	 	  (R0/R(PT))**(lk+1)*(-1.d0*PBAR1(K))
     +          *dble(-2-lk)*DCOS(DBLE(M)*LON(PT))
          AP(3,IS(M+1)+L+K-1)=1.d0/R(PT)**2*GM/R0*
CC     +          (R0/R(PT))**(lk+1)*PBAR1(K)*sin(LAT(PT))
     +		(R0/R(PT))**(lk+1)*(-1.d0*PBAR1(K))
     +          *dble(-2-lk)*DSIN(DBLE(M)*LON(PT))
          AP(4,IS(M+1)+K-1) =AP(2,IS(M+1)+K-1)
          AP(4,IS(M+1)+L+K-1)=AP(2,IS(M+1)+L+K-1)
          AP(5,IS(M+1)+K-1) =-1.d0/R(PT)**2*GM/R0*(R0/R(PT))**(lk+1)*
CC     +          (dble(lk+1)*PBAR(K)+cos(LAT(PT))*PBAR1(K)+
CC  use the next line,dcl-2010-2-13
     +	  (dble(lk+1)*PBAR(K)-Dcos(LAT(PT))*PBAR1(K)/Dsin(LAT(PT))+
     +    dble(m**2)*PBAR(K)/Dsin(LAT(PT))**2)*DCOS(DBLE(M)*LON(PT))
          AP(5,IS(M+1)+L+K-1)=-1.d0/R(PT)**2*GM/R0*(R0/R(PT))**(lk+1)*
CC     +          (dble(lk+1)*PBAR(K)+cos(LAT(PT))*PBAR1(K)+
CC  use the next line,dcl-2010-2-13
     +    (dble(lk+1)*PBAR(K)-Dcos(LAT(PT))*PBAR1(K)/Dsin(LAT(PT))+
     +    dble(m**2)*PBAR(K)/Dsin(LAT(PT))**2)*DSIN(DBLE(M)*LON(PT))
          AP(6,IS(M+1)+K-1) = -1.d0/(R(PT)**2*Dsin(LAT(PT)))*GM/R0*
     +          (R0/R(PT))**(lk+1)
     +          *PBAR(K)*dble(m)*dble(2+lk)*DSIN(DBLE(M)*LON(PT))
          AP(6,IS(M+1)+L+K-1)=+1.d0/(R(PT)**2*Dsin(LAT(PT)))*GM/R0*
     +          (R0/R(PT))**(lk+1)
     +          *PBAR(K)*dble(m)*dble(2+lk)*DCOS(DBLE(M)*LON(PT))
          AP(7,IS(M+1)+K-1) =AP(3,IS(M+1)+K-1)
          AP(7,IS(M+1)+L+K-1)=AP(3,IS(M+1)+L+K-1)
          AP(8,IS(M+1)+K-1) =AP(6,IS(M+1)+K-1)
          AP(8,IS(M+1)+L+K-1)=AP(6,IS(M+1)+L+K-1)
          AP(9,IS(M+1)+K-1) =GM/R0**3*(R0/R(PT))**(lk+3)
     +          *dble((lk+1)*(lk+2))*PBAR(K)*DCOS(DBLE(M)*LON(PT))
          AP(9,IS(M+1)+L+K-1)=GM/R0**3*(R0/R(PT))**(lk+3)
     +          *dble((lk+1)*(lk+2))*PBAR(K)*DSIN(DBLE(M)*LON(PT))
CC	  AP(1,IS(M+1)+K-1) = -AP(5,IS(M+1)+K-1)-AP(9,IS(M+1)+K-1)
CC	  AP(1,IS(M+1)+L+K-1)=-AP(5,IS(M+1)+L+K-1)-AP(9,IS(M+1)+L+K-1)
        ENDDO
        DEALLOCATE(PBAR,PBAR1,PBAR2)
        ENDDO
