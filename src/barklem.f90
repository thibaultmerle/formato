!
! From P. Barklem Version 2.0 http://www.astro.uu.se/~barklem/
!
Subroutine RETCROSS(NSlow,NSupp,Llow,Lupp,CROSS,ALPHA,IFAIL)
    !**********************************************************************
    !  Returns the cross section at 10000 m/s from tables for neutrals etc
    !   returns IFAIL=1 if not on the tables or some other error, else 0
    !  Paul Barklem Dec 1997
    !**********************************************************************
    !
    ! NSupp,NSlow = upper and lower state effective principal quantum numbers
    ! Lupp,Llow = angular momentum quantum numbers
    ! CROSS = hydrogen broadening cross-section for v=10000m/s
    ! ALPHA = velocity parameter
    !
    implicit none
    Real, intent(in) :: NSupp, NSlow
    Integer, intent(in) :: Lupp,Llow
    Real, intent(out):: CROSS,ALPHA
    Integer, intent(out) :: IFAIL
    Integer :: CHECK, TABLE
    !
    CHECK=ABS(Lupp-Llow)
    !
    if (CHECK.ne.1) then
        IFAIL=1
        !TM print *,'#######Forbidden transition!'
        return
    end if
    !
    TABLE=Llow+Lupp
    !
    if (TABLE.eq.1) call SP(NSlow,NSupp,Llow,Lupp,CROSS,ALPHA,IFAIL)
    if (TABLE.eq.3) call PD(NSlow,NSupp,Llow,Lupp,CROSS,ALPHA,IFAIL)
    if (TABLE.eq.5) call DF(NSlow,NSupp,Llow,Lupp,CROSS,ALPHA,IFAIL)
    if ((TABLE.ne.1).and.(TABLE.ne.3).and.(TABLE.ne.5)) then
        IFAIL=1
        !TM print *,'#######Transition not covered by current data'
        return
    end if
    !
    return
End Subroutine RETCROSS
!
!
!
Subroutine SP(NSlow,NSupp,Llow,Lupp,CROSS,ALPHA,IFAIL)
    !**********************************************************************
    !  Returns the cross section at 10000 m/s from s-p table 
    !   returns IFAIL=1 if not on the tables, else 0
    !  Paul Barklem Dec 1997
    !**********************************************************************
    !
    implicit none
    real*4 NSupp,NSlow,CROSS,ALPHA,NSS,NSP, &
               CS(21,18),AL(21,18),CS2(21,18),AL2(21,18), &
               NSSG(21),NSPG(18)
    integer Lupp,Llow,IFAIL,I,J
    character*70 COMMENTS
    !
    if (Lupp.eq.1) then
        NSP=NSupp
        NSS=NSlow
    else
        NSP=NSlow
        NSS=NSupp
    end if  
    !
    if ((NSS.gt.3.).or.(NSS.lt.1.).or.(NSP.gt.3.).or.(NSP.lt.1.3)) then
        IFAIL=1
        !TM print *,'#######Outside range of tables! (s-p)'
        return 
    end if
    !
    !  read in the table data
    !  3 lines of comments preceed
    !
    open(15,file='/home/tmerle/development/formato2/ad/bp/hb/sp.dat',form='formatted', status='old')
    !
    do I=1,3
        read(15,'(A70)') COMMENTS
    !   print *,COMMENTS
    end do   
    !
    do I=1,21
        read(15,*) (CS(I,J),J=1,18)
    end do
    !
    read(15,'(A70)') COMMENTS  ! a single comment line between data
    !
    do I=1,21
        read(15,*) (AL(I,J),J=1,18)
    end do
    !
    do I=1,21
        NSSG(I)=1.0+(I-1)*0.1
    end do
    !
    do I=1,18
        NSPG(I)=1.3+(I-1)*0.1
    end do
    !
    !   setup second derivative table for spline
    !
    call SPLIE2(NSSG,NSPG,CS,21,18,CS2)
    call SPLIE2(NSSG,NSPG,AL,21,18,AL2) 
    !
    !   run bicubic spline interpolation
    !
    call SPLIN2(NSSG,NSPG,CS,CS2,21,18,NSS,NSP,CROSS)
    call SPLIN2(NSSG,NSPG,AL,AL2,21,18,NSS,NSP,ALPHA)
    !
    close(15)
    return
    !
End Subroutine SP
!
!
!
Subroutine PD(NSlow,NSupp,Llow,Lupp,CROSS,ALPHA,IFAIL)
    !**********************************************************************
    !  Returns the cross section at 10000 m/s from p-d table 
    !   returns IFAIL=1 if not on the tables, else 0
    !  Paul Barklem Dec 1997
    !**********************************************************************
    !
    implicit none
    real*4 NSupp,NSlow,CROSS,ALPHA,NSD,NSP, &
              CS(18,18),AL(18,18),CS2(18,18),AL2(18,18), &
              NSPG(18),NSDG(18)
    integer Lupp,Llow,IFAIL,I,J
    character*70 COMMENTS
    !
    if (Lupp.eq.2) then
        NSD=NSupp
        NSP=NSlow
    else
        NSD=NSlow
        NSP=NSupp
    end if
    !
    if ((NSP.gt.3.).or.(NSP.lt.1.3).or.(NSD.gt.4.).or.(NSD.lt.2.3)) then
        IFAIL=1
        !TM print *,'#######Outside range of tables! (p-d)'
        return 
    end if
    !
    !  read in the table data
    !  3 lines of comments preceed
    !
    open(15,file='/home/tmerle/development/formato2/ad/bp/hb/pd.dat',form='formatted', &
        status='old')
    !
    do I=1,3
        read(15,'(A70)') COMMENTS
    !   print *,COMMENTS
    end do       
    !
    do I=1,18
        read(15,*) (CS(I,J),J=1,18)
    end do
    !
    read(15,'(A70)')COMMENTS  ! a single comment line between data
    !
    do I=1,18
        read(15,*) (AL(I,J),J=1,18)
    end do
    !
    !  make arrays of the nstar values
    !
    do I=1,18
        NSPG(I)=1.3+(I-1)*0.1
        NSDG(I)=2.3+(I-1)*0.1
    end do
    !
    !   setup second derivative table for spline
    !
    call SPLIE2(NSPG,NSDG,CS,18,18,CS2)
    call SPLIE2(NSPG,NSDG,AL,18,18,AL2) 
    !
    !   run bicubic spline interpolation
    !
    call SPLIN2(NSPG,NSDG,CS,CS2,18,18,NSP,NSD,CROSS)
    call SPLIN2(NSPG,NSDG,AL,AL2,18,18,NSP,NSD,ALPHA)
    !
    !
    close(15)
    return
    !
End Subroutine PD
!
!
!
Subroutine DF(NSlow,NSupp,Llow,Lupp,CROSS,ALPHA,IFAIL)
    !**********************************************************************
    !  Returns the cross section at 10000 m/s from d-f table 
    !   returns IFAIL=1 if not on the tables, else 0
    !  Paul Barklem Dec 1997
    !**********************************************************************
    !
    implicit none
    real*4 NSupp,NSlow,CROSS,ALPHA,NSD,NSF, &
               CS(18,18),AL(18,18),CS2(18,18),AL2(18,18), &
               NSDG(18),NSFG(18)
    integer Lupp,Llow,IFAIL,I,J
    character*70 COMMENTS
    !
    if (Lupp.eq.3) then
        NSF=NSupp
        NSD=NSlow
    else
        NSF=NSlow
        NSD=NSupp
    end if
    !
    if ((NSD.gt.4.).or.(NSD.lt.2.3).or.(NSF.gt.5.).or.(NSF.lt.3.3)) then
        IFAIL=1
        !TM print *,'#######Outside range of tables! (d-f)'
        return 
    end if
    !
    !  read in the table data
    !  3 lines of comments preceed
    !
    open(15,file='/home/tmerle/development/formato2/ad/bp/hb/df.dat',form='formatted', status='old')
    !
    do I=1,3
        read(15,'(A70)') COMMENTS
    !   print *,COMMENTS
    end do       
    !
    do I=1,18
        read(15,*) (CS(I,J),J=1,18)
        !print *,(CS(I,J),J=1,18)
    end do
    !
    read(15,'(A70)')COMMENTS  ! a single comment line between data
    !
    do I=1,18
        read(15,*) (AL(I,J),J=1,18)
        !print *,(AL(I,J),J=1,18)
    end do
    !
    !  make arrays of the nstar values
    !
    do I=1,18
        NSDG(I)=2.3+(I-1)*0.1
        NSFG(I)=3.3+(I-1)*0.1
    end do
    !
    !   setup second derivative table for spline
    !
    call SPLIE2(NSDG,NSFG,CS,18,18,CS2)
    call SPLIE2(NSDG,NSFG,AL,18,18,AL2) 
    !
    !   run bicubic spline interpolation
    !
    call SPLIN2(NSDG,NSFG,CS,CS2,18,18,NSD,NSF,CROSS)
    call SPLIN2(NSDG,NSFG,AL,AL2,18,18,NSD,NSF,ALPHA)
    !
    !
    close(15)
    return
    !
End Subroutine DF
!
!**********************************************************************
!  The following are Numerical Recipes routines
!   W. Press, S. Teukolsky, W. Vetterling and B. Flannery
!  " Numerical Recipes in Fortran: The Art of Scientific 
!     Computing"   Second Edition, Cambridge Univ. Press
!**********************************************************************
!
!----------------------------------------------------------------------
!
      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)
      PARAMETER (NN=100)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN)
      DO 13 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
        CALL SPLINE(X2A,YTMP,N,1.E30,1.E30,Y2TMP)
        DO 12 K=1,N
          Y2A(J,K)=Y2TMP(K)
12      CONTINUE
13    CONTINUE
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
!
      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
      PARAMETER (NN=100)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN), &
                YYTMP(NN)
      DO 12 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
          Y2TMP(K)=Y2A(J,K)
11      CONTINUE
        CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J))
12    CONTINUE
      CALL SPLINE(X1A,YYTMP,M,1.E30,1.E30,Y2TMP)
      CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y)
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
!
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
           /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
!
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
!     TM remplace PAUSE par STOP
      IF (H.EQ.0.) STOP 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+ &
           ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
