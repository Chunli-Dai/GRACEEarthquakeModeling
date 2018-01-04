        M=0
        L=NMAX+1
        ALLOCATE(PBAR(L),PBAR1(L),PBAR2(L))
        CALL LGDR2(DCOS(LAT(PT)),NMAX,M,PBAR,PBAR1,PBAR2)
        DO K=1,L
CC        PNMC(IS(1)+K-1) = (R0/R(PT))**(M+K) * (PBAR(K))
          lk=k+M-1
          AP(3,IS(1)+K-1) =1.D0/R(PT)*GM/R0*(R0/R(PT))**(lk+1)*
     +		(dble(lk+1)*PBAR(K))
          AP(1,IS(1)+K-1) =-1.d0/R(PT)*GM/R0*(R0/R(PT))**(lk+1)
     +		*(PBAR1(K))
          AP(2,IS(1)+K-1) =0.d0
        ENDDO
        DEALLOCATE(PBAR,PBAR1,PBAR2)

        DO M=1,NMAX
!       DO M=1,3  !to test
        L=NMAX-M+1
        ALLOCATE(PBAR(L),PBAR1(L),PBAR2(L))
        CALL LGDR2(DCOS(LAT(PT)),NMAX,M,PBAR,PBAR1,PBAR2)
        DO K=1,L
CC      PNMC(IS(M+1)+K-1  ) = (R0/R(PT))**(M+K) * (PBAR(K)) * COS(M*LON(PT))
CC      PNMC(IS(M+1)+L+K-1) = (R0/R(PT))**(M+K) * (PBAR(K)) * SIN(M*LON(PT))
          lk=k+M-1      ! The number of l
!-Vr
          AP(3,IS(M+1)+K-1) =1.D0/R(PT)*GM/R0*(R0/R(PT))**(lk+1)
     +		*(dble(lk+1)*PBAR(K))*DCOS(dble(M)*LON(PT))
          AP(3,IS(M+1)+L+K-1)=1.D0/R(PT)*GM/R0*(R0/R(PT))**(lk+1)
     +		*(dble(lk+1)*PBAR(K))*DSIN(DBLE(M)*LON(PT))
!-1/r*Vsita
          AP(1,IS(M+1)+K-1) = 1.d0/R(PT)*GM/R0*
     +	 	  (R0/R(PT))**(lk+1)*(-1.d0*PBAR1(K))
     +          *DCOS(DBLE(M)*LON(PT))
          AP(1,IS(M+1)+L+K-1)=1.d0/R(PT)*GM/R0*
     +		(R0/R(PT))**(lk+1)*(-1.d0*PBAR1(K))
     +          *DSIN(DBLE(M)*LON(PT))
!1/r/sin(colat)*Vlamda
          AP(2,IS(M+1)+K-1) = -1.d0/(R(PT)*Dsin(LAT(PT)))*GM/R0*
     +          (R0/R(PT))**(lk+1)
     +          *PBAR(K)*dble(m)*DSIN(DBLE(M)*LON(PT))
          AP(2,IS(M+1)+L+K-1)=+1.d0/(R(PT)*Dsin(LAT(PT)))*GM/R0*
     +          (R0/R(PT))**(lk+1)
     +          *PBAR(K)*dble(m)*DCOS(DBLE(M)*LON(PT))
        ENDDO
        DEALLOCATE(PBAR,PBAR1,PBAR2)
        ENDDO
