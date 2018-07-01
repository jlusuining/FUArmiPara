
program FindOutbursts
  implicit none

  integer, parameter :: DP=8
  integer, parameter :: nar=5
  integer :: iost, i, nlt, ngt, noutb

  integer :: nmid

  real(DP) :: t(1:nar), Mdot(1:nar), Mdisk(1:nar), Mstar(1:nar)

  real(DP) :: toutburst(1000)

  real(DP) :: Mdotcrit=1.0E-4_DP


  if (mod(nar,2)==1) then
    nmid=nar/2+1
  else
    nmid=nar/2
  end if
  write(*,*) 'Nmid=', nmid

  toutburst=0.0_DP

  open(1,file='m.txt',status='old')
  read(1,*) !'#line'
  do i=1, nar
    read(1,*, IOSTAT=iost) t(i), Mdot(i), Mdisk(i), Mstar(i)
    if (iost/=0) exit
  end do

  noutb=0

  do
    call backvalues(t)
    call backvalues(Mdot)
    call backvalues(Mdisk)
    call backvalues(Mstar)

    read(1,*, IOSTAT=iost) t(nar), Mdot(nar), Mdisk(nar), Mstar(nar)
    if (iost/=0) exit

    !do i=1, 5
    !  write(*,*) i, t(i), Mdot(i), Mdisk(i), Mstar(i)
    !enddo
    !read(*,*)

    nlt=0
    do i=1, nmid-1
      if (Mdot(i)<Mdotcrit) then
        nlt=nlt+1
      end if
    end do

    ngt=0
    do i=nmid, nar
      if (Mdot(i)>Mdotcrit) then
        ngt=ngt+1
      end if
    end do

    if (nlt==nmid-1) then
      !left all < Mdotcrit
      if (ngt==nar-nmid+1) then
        !right all > Mdotcrit
        !found 1
        noutb=noutb+1
        toutburst(noutb)=t(nmid)
        write(*,*) noutb, toutburst(noutb)
      end if

    end if


  end do


  open(2,file='toutburst.txt')

  do i=1, noutb
    write(2,*) toutburst(i)
  end do

  do i=1, noutb-1
    call subtitute(i, toutburst(i)-200.0_DP, toutburst(i+1))
  end do


contains
  subroutine backvalues(arr)
    !arr(nar)=arr(nar-1), rr(4)=arr(3), rr(3)=arr(2), rr(2)=arr(1), rr(1)left
    implicit none
    real(DP), intent(inout) :: arr(1:nar)
    integer :: i

    do i=2, nar
      arr(i-1)=arr(i)
    end do

  end subroutine backvalues

  subroutine subtitute(io, to1, to2)
    implicit none
    real(DP), intent(in) :: to1, to2
    integer, intent(in) :: io
    integer :: ios
    character(100) :: st, st2, filename
    character(10) :: sto1, sto2
    character(7) :: sio

    if (io<10) then
      write(sio, '(A2, I1)') '00', io
    else if (io<100) then
      write(sio, '(A1, I2)') '0', io
    else if (io<1000) then
      write(sio, '(I3)') io
    else
      write(*,*) 'Can not handel filenumber >= 1000'
      stop
    end if
    sio=adjustl(trim(sio))//'full'

    write(sto1, '(F10.1)') to1
    sto1=adjustl(sto1)
    write(sto2, '(F10.1)') to2
    sto2=adjustl(sto2)

    write(filename,'(A4, A3, A4)') 'mdot', sio, '.gnu'
    filename=adjustl(trim(filename))
    open(4,file=filename)
    open(3,file='mdot_detail_samp.gnu')
    do
      read(3,'(A)',iostat=ios) st
      st2=Replace_Text(st, '[filenumber]', sio)
      st=Replace_Text(st2, '[tout1]', sto1)
      st2=Replace_Text(st, '[tout2]', sto2)

      write(4,'(A)') st2(1:len_trim(st2))
      if (ios/=0) exit
      !write(*,*) st
    end do
    close(3)
    close(4)
  end subroutine subtitute

  FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
    CHARACTER(*)        :: s,text,rep
    CHARACTER(LEN(s)+100) :: outs     ! provide outs with extra 100 char len
    INTEGER             :: i, nt, nr

    outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
    DO
      i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
      outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
    END DO
  END FUNCTION Replace_Text
end program

