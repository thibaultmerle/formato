Subroutine format_px(flag, ifile, ofile_idx, ofile_px)
    !    
    Implicit None
    !
    Character(len=*), intent(in) :: flag
    Character(len=*), intent(in) :: ifile
    Character(len=*), intent(in) :: ofile_idx
    Character(len=*), intent(in) :: ofile_px
    !
    Integer :: i, ios, npx_tot = 0, npx_dtl
    Integer :: mult, L, parity, pos, counter = 0
    Character(len=256) :: line
    Character(len=3) :: cterm
    Real, dimension(:), allocatable :: e_ryd, x_Mb 
    Real :: dumb
    Integer :: dum1, dum2, dum3
    !
    Open(unit=10, file=ifile, status='old', action='read', position='rewind', iostat=ios)
    IF (ios /= 0) STOP "Pb with the opening of the photoionization file."
    !
    Open(unit=20, file=ofile_idx, status='replace', action='write', form='formatted', iostat=ios)
    IF (ios /= 0) STOP "Pb with the opening of the index file."
    !
    Open(unit=30, file=ofile_px, status='replace', action='write', access='direct', recl=8, iostat=ios)
    IF (ios /= 0) STOP "Pb with the opening of the photoionization binary file."
    !
    Do 
        Read(10, '(A)', iostat=ios) line
        IF (ios /= 0) STOP "Pb with the reading of comment lines."
        IF (Scan(line, '#') == 0) EXIT
    End Do
    !
    Backspace(10)
    !
    Do
        counter = counter + npx_tot + 1
        If (flag == 'norad') Then
            Read(10, *, iostat=ios) mult, L, parity, pos
            IF (ios > 0) STOP "Pb with the reading of mult, L, parity, pos"
            IF (ios < 0) EXIT
            !        
            Read(10, *, iostat=ios) npx_dtl, npx_tot
            If (ios /= 0) Then
                Backspace(10)
                Read(10, *, iostat=ios) dumb, npx_tot
                Write(20, *) counter, dumb, npx_tot+1, mult, L, parity, pos
            Else
                Write(20, *) counter, npx_dtl, npx_tot+1, mult, L, parity, pos
            End If
        Else If (flag == 'topbase') Then
            Read(10, *, iostat=ios) dum1, dum2, dum3, cterm, pos, dumb, npx_tot 
            IF (ios > 0) STOP "Pb with the reading of dumb, dumb, dumb, mult, pos, dumb, npx_tot"
            IF (ios < 0) EXIT
            
            !Write(*, *) ios, dum1, dum2, dum3, cterm, pos, dumb, npx_tot
            npx_tot = npx_tot -1
            Read(cterm(1:1), '(I1)') mult
            Read(cterm(2:2), '(I1)') L 
            Read(cterm(3:3), '(I1)') parity
            Write(20, *) counter, dumb, npx_tot+1, mult, L, parity, pos
        Else
            STOP "Only flags allowed 'norad' or 'topbase'"
        End If
        !
        ! Empty table
        If (npx_tot+1 == 0) Then
            CYCLE
        End If
        !
        Allocate(e_ryd(npx_tot+1), x_Mb(npx_tot+1), stat=ios)
        IF (ios /= 0) STOP "Allocation memory pb."
        e_ryd = 0.
        x_Mb = 0.
        !
        Read(10, *, iostat=ios) (e_ryd(i), x_Mb(i), i=1, npx_tot+1)
        IF (ios /= 0) STOP "Pb with the reading of a table."
        !
        Do i=1, npx_tot+1
            Write(30, rec=counter+i-1) e_ryd(i), x_Mb(i)
        End Do
        !
        Deallocate(e_ryd)
        Deallocate(x_Mb)
    End Do
    !
    Close(10)
    Close(20)
    Close(30)
    !
    Write(*, *) "Read file ok."
    Write(*, *) "File ", ofile_idx, " created."
    Write(*, *) "File ", ofile_px, " created."
    !
End Subroutine format_px
!
Subroutine extract_px(ifile, idx, n, e_ryd, x_Mb)
    !
    Implicit None
    !
    Character(len=*), intent(in) :: ifile
    Integer, intent(in) :: idx
    Integer, intent(in) :: n
    Real, intent(out), dimension(n) :: e_ryd
    Real, intent(out), dimension(n) :: x_Mb
    Integer :: i, ios
    !
    e_ryd = 0.; x_Mb = 0.
    !
    !Write(*, *) "idx = ", idx, "n = ", n
    !
    Open(unit=10, file=ifile, status='old', action='read', access='direct',recl=8, iostat=ios)
    IF (ios /= 0) STOP "Pb with the opening of the binary photoionization file." 
    !
    Do i=idx, idx+n-1
        Read(10, rec=i) e_ryd(i+1-idx), x_Mb(i+1-idx) 
    End Do
    !
    Close(10)
    !
End Subroutine extract_px
!
Program test
    !
    Implicit None
    !
    !Character(len=*), parameter :: flag = 'norad'
    !Character(len=*), parameter :: ifile = '../ad/rt/bf/norad_fei_px.dat' 
    !Character(len=*), parameter :: ofile_idx = '../ad/rt/bf/norad_fei_idx.dat' 
    !Character(len=*), parameter :: ofile_px = '../ad/rt/bf/norad_fei_px.bin' 
    !
    Character(len=*), parameter :: flag = 'topbase'
    Character(len=*), parameter :: ifile = '../ad/rt/bf/topbase_mgi_px.dat' 
    Character(len=*), parameter :: ofile_idx = '../ad/rt/bf/topbase_mgi_idx.dat' 
    Character(len=*), parameter :: ofile_px = '../ad/rt/bf/topbase_mgi_px.bin'   
    !
    Integer :: i, idx, n 
    Real, dimension(10000) :: e, x
    !
    Call format_px(flag, ifile, ofile_idx, ofile_px)
    !
    !idx = 1; n = 2514
    idx = 1; n = 743
    Call extract_px(ofile_px, idx, n, e, x)
    Write(*, '(ES12.6, ES10.3)') (e(i), x(i), i=idx, 10)
    Write(*, *) '...'
    Write(*, '(ES12.6, ES10.3)') (e(i), x(i), i=n-2, n)
    !
    !Call extract_px(ofile_px, 2515, 2601, e, x)
    !Call extract_px(ofile_px, 5116, 2628, e, x)   
    !
End Program test