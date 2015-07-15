!
! From Janicki, CoPhyC (1990)
!
Subroutine facln(N, result)
    !
    ! From Numerical recipes
    !
    Implicit None
    Integer, Intent(in) :: N
    Double Precision, intent(out) :: result 
    Double Precision, Parameter :: stp = 2.506628227465D0
    Double Precision, Dimension(6), Parameter :: &
    cof = (/76.18009173D0, -86.50532033D0, 24.01409822D0, -1.231739516D0, 0.120858003D-2, -0.536382D-5/)
    Double Precision :: x, tmp, ser
    Integer :: I
    !
    x = Dble(N)
    !
    tmp = (x+0.5D0) * log(x+5.5D0) - (x+5.5D0)
    !
    ser = 1.D0
    !
    Do I = 1, 6
        x = x + 1.D0
        ser = ser + cof(I)/x
    End Do
    !
    result = tmp + Log(stp*ser)
    !
    Return
    !
End Subroutine facln
!
Subroutine gbf_ave(N, Z, hv, gbf_all, gbf_av)
    !
    Implicit None
    !
    Integer, intent(in) :: N, Z
    Double Precision, intent(in) :: hv
    Double Precision, dimension(N), intent(out) :: gbf_all
    Double Precision, intent(out) :: gbf_av
    Double Precision :: gbf
    Integer :: L
    !
    gbf_av = 0.D0
    !
    Do L = 0, N-1
        Call gaunt_bf(N, L, Z, hv, gbf)
        gbf_all(L+1) = gbf
        gbf_av = gbf_av + (2*L+1)*gbf
    End Do
    !
    gbf_av = gbf_av/Dble(N**2)
    !
    Return
    !
End Subroutine gbf_ave
!
Subroutine gaunt_bf(N, L, Z, hv, gbf)
    !
    Implicit None
    !
    Integer, intent(in) :: N, L, Z
    Double Precision, intent(in) :: hv
    Double Precision, intent(out) :: gbf
    Double Precision, Parameter :: Ryd = 13.6057D0
    Double Precision, Parameter :: Pi = 3.141592653D0
    Double Precision :: Ee, rho, eta, R12, EL2, A1, A2, A3, A4, sum
    Double Precision :: facln1, facln2, facln3, facln4
    Double Precision :: gl1, gl2, gl3, gl4
    Integer :: I
    !
    ! Original form but cannot reproduce properly Fig.3 of Janicki 1990
    !Ee = hv - Z**2.D0 * Ryd / Dble(N**2)
    !
    ! This works
    Ee = hv
    !
    !Write(*, *) "Ee:", Ee
    !
    gbf = 0.D0
    !
    IF (Ee <= 1.D-18) RETURN
    !
    eta = sqrt(Z**2.D0 * Ryd / Ee)
    rho = eta/Dble(N)
    !
    !Write(*,*) "eta, rho:", eta, rho
    !
    sum = 0.D0
    !
    Do I = 1, L+1
        sum = sum + Log(I**2.D0 + eta**2.D0)
    End Do
    !
    R12 = 1.D0 + rho**2.D0
    EL2 = eta**2.D0 + (L+1)**2.D0
    !Write(*, *) "R12, EL2:", R12, EL2
    !
    Call facln(N+L, facln1)
    Call facln(2*L+1, facln2)
    Call facln(2*L+2, facln3)
    Call facln(N-L-1, facln4)
    !
    !Write(*, *) "facln:", facln1, facln2, facln3, facln4
    !
    A1 = Exp(-2.D0*eta*(Pi/2.D0 - Atan(rho)) + &
             0.5D0*(facln1 - facln2 - facln3 - facln4 + sum + L*Log(16.D0)) + &
    &        (L-1.D0)*Log(rho) - (N-2.D0)*Log(R12))

    A2 = (0.3401D0*N/EL2) / (1.D0 - Exp(-2.D0*Pi*eta))
    A3 = 4.D0*L**3.D0 * (L+1.D0) * (2*L+1.D0) / (L**2.D0+eta**2.D0)
    A4 = 64.D0 * (L+1.D0)**2.D0 * (rho*eta/R12)**2.D0 / ((2.D0*L+1.D0)*EL2)
    !
    !Write(*, *) "A:", A1, A2, A3, A4
    !
    Call gl_bf(L, L+1-N, eta, rho, gl1)
    Call gl_bf(L, L-1-N, eta, rho, gl2)
    Call gl_bf(L+1, L+1-N, eta, rho, gl3)
    Call gl_bf(L+1, L-N, eta, rho, gl4)
    !
    !Write(*, *) "gl:", gl1, gl2, gl3, gl4
    !
    gbf = A2*(A3*(A1*gl1-A1*gl2/R12**2.D0)**2.D0 + A4*((L+1-N)*A1*gl3 + A1*gl4*(L+1+N)/R12)**2.D0)
    !
    Return
    !
End Subroutine
!
Subroutine gl_bf(L, M_minus, eta, rho, gl)
    !
    Implicit None
    !
    Integer, intent(in) :: L, M_minus
    Double Precision, intent(in) :: eta, rho
    Double Precision, intent(out) :: gl
    Double Precision, dimension(210) :: b
    Integer :: I, M
    !
    gl = 0.D0
    !
    IF (L == 0) RETURN
    !
    M = -M_minus
    b(1) = 1.0D0
    b(2) = 2 * M * eta / Dble(L)
    !
    Do I=2, 2*M
        b(I+1) = -(4*eta*(I-1-M)*b(I)+(2*M+2-I)*(2*M+2*L+1-I)*B(I-1)) / (I*(I+2*L-1))
        !Write(*, *)"b:", b(I+1)
    End Do
    !
    Do I = 0, 2*M
        gl = gl + b(I+1)*rho**I
    End Do
    !
    Return
    !
End Subroutine gl_bf
