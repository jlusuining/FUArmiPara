module ConstantMod
  implicit none
  !variable type 8 for double precision
  integer, parameter :: DP=8
  !constants
  real(DP), parameter :: Pi=3.14159265358979_DP
  real(DP), parameter :: Grav=6.672E-8_DP    !cgs
  real(DP), parameter :: Rgas=8.31441E7_DP   !cgs
  real(DP), parameter :: Stfb=5.67032E-5_DP  !cgs
  real(DP), parameter :: Msun=1.9891E33_DP   !cgs
  real(DP), parameter :: syear=31556952_DP   !cgs
  real(DP), parameter :: AU=1.49597870E13_DP !cgs
  real(DP), parameter :: Lsun=3.9E33_DP      !cgs erg s-1
  real(DP), parameter :: Rsun=6.955E10_DP    !cgs cm
end module

module input_parameter_mod
  use ConstantMod
  implicit none
  !variable type 8 for double precision
  !integer, parameter :: DP=8

  real(DP) :: Rin_input
  real(DP) :: Rout_input
  real(DP) :: tend_input

  real(DP) :: Mdotinfall_input
  real(DP) :: Radd_input
  real(DP) :: dRadd_input

  real(DP) :: alpha_MRI
  real(DP) :: alpha_GI0
  real(DP) :: alpha_Dead

  real(DP) :: Tcrit_input
  real(DP) :: Sigcrit_input
  real(DP) :: Qcirt_input

  real(DP) :: Mstar_input

  real(DP) :: tout_input
  real(DP) :: percent_Mdot_input
  real(DP) :: outMdot

  !gives the average Mdot(r) from radial zone iMds to iMde
  !Mdaver = \Sum (Mdr(iMds), Mdr(iMde)) / (iMde-iMds+1)
  integer :: iMds, iMde

contains
  subroutine get_parameters()
    implicit none
    open(1,file='inpara_alpdead.txt',status='old')
    read(1,*) Rin_input
    read(1,*) Rout_input
    read(1,*) tend_input


    read(1,*) Mdotinfall_input
    read(1,*) Radd_input
    read(1,*) dRadd_input

    read(1,*) alpha_MRI
    read(1,*) alpha_GI0
    read(1,*) alpha_Dead

    read(1,*) Tcrit_input
    read(1,*) Sigcrit_input
    read(1,*) Qcirt_input

    read(1,*) Mstar_input

    read(1,*) tout_input
    !read(1,*) percent_Mdot_input
    read(1,*) outMdot
    read(1,*) iMds, iMde

    close(1)

    write(*,*) Rin_input, 'Rin cm'
    write(*,*) Rout_input, 'Rout cm'
    write(*,*) tend_input, 'tEnd year'


    write(*,*) Mdotinfall_input, 'Mdotinfall Msun/year'
    write(*,*) Radd_input, 'Radd AU'
    write(*,*) dRadd_input, 'dRadd AU'

    write(*,*) alpha_MRI, 'alpha_MRI'
    write(*,*) alpha_GI0, 'alpha_GI0'
    write(*,*) alpha_Dead, 'alpha_Dead'

    write(*,*) Tcrit_input, 'Tcrit K'
    write(*,*) Sigcrit_input, 'Sig_crit g cm-2 (default 200)'
    write(*,*) Qcirt_input, 'Qcrit default 2'

    write(*,*) Mstar_input, 'Mstar initail (Msun)'

    write(*,*) tout_input, 'tout Sig Tem every ... year'
    !write(*,*) percent_Mdot_input, 'percent_Mdot'
    write(*,*) outMdot, 'out put Mdot every ... year'
    write(*,*) iMds, iMde, 'average the Mdot(r) from radial bin iMds to iMde'
  end subroutine

end module


MODULE tridag_mod
  !module for sovling tridag
  !from Numerical recipes
  implicit none
contains

  SUBROUTINE tridag_ser(a,b,c,r,u)
    IMPLICIT NONE
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(8), DIMENSION(:), INTENT(OUT) :: u
    REAL(8), DIMENSION(size(b)) :: gam
    INTEGER(4) :: n,j
    REAL(8) :: bet
    n=assert_eqn((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
    bet=b(1)
    if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
    u(1)=r(1)/bet
    do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j-1)*gam(j)
      if (bet == 0.0) &
        call nrerror('tridag_ser: Error at code stage 2')
      u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridag_ser

  RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
    IMPLICIT NONE
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(8), DIMENSION(:), INTENT(OUT) :: u
    INTEGER(4), PARAMETER :: NPAR_TRIDAG=4
    INTEGER(4) :: n,n2,nm,nx
    REAL(8), DIMENSION(size(b)/2) :: y,q,piva
    REAL(8), DIMENSION(size(b)/2-1) :: x,z
    REAL(8), DIMENSION(size(a)/2) :: pivc
    n=assert_eqn((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
    if (n < NPAR_TRIDAG) then
      call tridag_ser(a,b,c,r,u)
    else
      if (maxval(abs(b(1:n))) == 0.0) &
        call nrerror('tridag_par: possible singular matrix')
      n2=size(y)
      nm=size(pivc)
      nx=size(x)
      piva = a(1:n-1:2)/b(1:n-1:2)
      pivc = c(2:n-1:2)/b(3:n:2)
      y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
      q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
      if (nm < n2) then
        y(n2) = b(n)-piva(n2)*c(n-1)
        q(n2) = r(n)-piva(n2)*r(n-1)
      end if
      x = -piva(2:n2)*a(2:n-2:2)
      z = -pivc(1:nx)*c(3:n-1:2)
      call tridag_par(x,y,z,q,u(2:n:2))
      u(1) = (r(1)-c(1)*u(2))/b(1)
      u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
        -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
      if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
    end if
  END SUBROUTINE tridag_par

  !BL
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
      assert_eqn=nn(1)
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag:', &
        string
      STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn
  !BL
  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror
  !BL
END MODULE



module CloudCoreMod
  use ConstantMod
  use input_parameter_mod
  implicit none
  real(DP) :: Temp_Core
  real(DP) :: Radius_Core
  real(DP) :: W_Core
  real(DP) :: Mass_Core
  real(DP) :: Mu_Core
  real(DP) :: tlife_Core

contains
  function SoundSpeed2_Core()
    implicit none
    real(DP) :: SoundSpeed2_Core
    !implicit function of Temp_Core
    SoundSpeed2_Core=Rgas*Temp_Core/Mu_Core
  end function SoundSpeed2_Core

  function Mdot_Core()
    implicit none
    real(DP) :: Mdot_Core
    real(DP) :: a3
    a3=SoundSpeed2_Core()**1.5_DP
    Mdot_Core=0.975_DP*a3/Grav
  end function Mdot_Core

  function Rd_Core(t)
    implicit none
    real(DP) :: Rd_Core
    real(DP), intent(in) :: t
    real(DP), save :: coe
    integer, save :: itimes
    if (itimes/=123456) then
      coe=SoundSpeed2_Core()**2.0_DP*W_Core*W_Core/16.0_DP/Grav/Mdot_Core()
      itimes=123456
    endif
    Rd_Core=coe*t*t*t
  end function Rd_Core

  function Srt_Core(R_disk, t, it)
    implicit none
    real(DP) :: Srt_Core
    real(DP), intent(in) :: R_disk, t
    integer, intent(in) :: it
    real(DP), save :: Mdotsave
    integer, save :: itsave
    real(DP), save :: Rdsave
    integer, save :: itimes
    if (it/=itsave.OR.itimes/=123456) then
      Rdsave=Rd_Core(t)
      Mdotsave=Mdot_Core()
      itimes=123456
      itsave=it
      !write(*,*) 'coe'
    endif
    if (R_disk<Rdsave.and.t<tlife_Core) then
      Srt_Core=Mdotsave/4.0_DP/Pi/R_disk/Rdsave/sqrt(1.0_DP-R_disk/Rdsave)
    else
      Srt_Core=0.0_DP
    endif
  end function Srt_Core

  function Srtj_Core(R_disk, t, it)
    implicit none
    real(DP) :: Srtj_Core
    real(DP) :: Srt_Core
    real(DP), intent(in) :: R_disk, t
    integer, intent(in) :: it
    real(DP), save :: Mdotsave
    integer, save :: itsave
    real(DP), save :: Rdsave
    integer, save :: itimes
    real(DP) :: x
    if (it/=itsave.OR.itimes/=123456) then
      Rdsave=Rd_Core(t)
      Mdotsave=Mdot_Core()
      itimes=123456
      itsave=it
      !write(*,*) 'coe'
    endif
    if (R_disk<Rdsave.and.t<tlife_Core) then
      x=R_disk/Rdsave
      Srt_Core=Mdotsave/4.0_DP/Pi/R_disk/Rdsave/sqrt(1.0_DP-R_disk/Rdsave)
      Srtj_Core=Srt_Core+Srt_Core*(2.0_DP-3.0_DP*sqrt(x)+x/(1.0_DP+sqrt(x)))
    else
      Srtj_Core=0.0_DP
    endif
  end function Srtj_Core


  !tech fun sub
  subroutine ReadCloudCorePara(parafname)
    implicit none
    !tech var
    character(*), intent(in) :: parafname
    character(20) :: st_temp
    integer :: iost
    !end tech var

    open(10, file=parafname, status='old')
    do
      read(10, *, iostat=iost) st_temp
      if (iost/=0) exit
      if (index(st_temp, 'Temp_Core')/=0) then
        read(10, *) Temp_Core
      elseif (index(st_temp, 'W_Core')/=0) then
        read(10, *) W_Core
      elseif (index(st_temp, 'Mass_Core')/=0) then
        read(10, *) Mass_Core
      elseif (index(st_temp, 'Mu_Core')/=0) then
        read(10, *) Mu_Core
      endif
    enddo
    close(10)
    Mass_Core=Mass_Core*Msun
    tlife_Core=Mass_Core/Mdot_Core()
    write(*,*) '--------Core Para--------'
    write(*,*) 'W_Core', W_Core
    write(*,*) 'Temp_Core', Temp_Core
    write(*,*) 'Mass_Core/Msun', Mass_Core/Msun
    write(*,*) 'Mu_Core', Mu_Core
    write(*,*) 'tlife_Core/syear', tlife_Core/syear
    write(*,*) '-------------------------'
  end subroutine ReadCloudCorePara

  function Mdot2Starfcore(t)
    !Mass falls onto the star from the core collapse
    !Eq 18 of Jin & Sui 2010
    implicit none
    real(DP), intent(in) :: t
    real(DP) :: Mdot2Starfcore
    if (t<tlife_Core) then
      Mdot2Starfcore=Mdot_Core()*(1.0_DP-sqrt(1.0_DP-Rin_input/Rd_Core(t)))
    else
      Mdot2Starfcore=0.0_DP
    end if
  end function
end module CloudCoreMod

module StellarMod
  use ConstantMod
  implicit none
  real(DP) :: Mass_Stellar
  real(DP) :: Radius_Stellar
  real(DP) :: Temp_Stellar
contains
  subroutine AddMass_Stellar(AddMass)
    implicit none
    real(DP), intent(in) :: AddMass
    Mass_Stellar=Mass_Stellar+AddMass
  end subroutine AddMass_Stellar

  !tech fun sub
  subroutine ReadStellarPara(parafname)
    implicit none
    !tech var
    character(*), intent(in) :: parafname
    character(20) :: st_temp
    integer :: iost
    !end tech var

    open(11, file=parafname, status='old')
    do
      read(11, *, iostat=iost) st_temp
      if (iost/=0) exit
      if (index(st_temp, 'Temp_Stellar')/=0) then
        read(11, *) Temp_Stellar
      elseif (index(st_temp, 'Radius_Stellar')/=0) then
        read(11, *) Radius_Stellar
      endif
    enddo
    close(11)
    Radius_Stellar=Radius_Stellar*AU
    write(*,*) '------Stellar Para--------'
    write(*,*) 'Temp_Stellar', Temp_Stellar
    write(*,*) 'Radius_Stellar/AU', Radius_Stellar/AU
    write(*,*) '-------------------------'
  end subroutine ReadStellarPara
end module StellarMod



module evolve_mod
  use input_parameter_mod
  !module for evolution of the disk
  use tridag_mod
  use CloudCoreMod
  implicit none

  !variable type 8 for double precision
  !integer, parameter :: DP=8
  !length of array (Martin=200)
  integer, parameter :: nar=1200 !Armitage

  !struct for disk
  type :: disk
    real(DP) :: t
    real(DP) :: R(0:nar+1), Sig(0:nar+1), Tc(0:nar+1)
    real(DP) :: Te(0:nar+1), Sigm(0:nar+1), Sigg(0:nar+1), QToom(0:nar+1)
    real(DP) :: Rin, Rout
    real(DP) :: Mu
    real(DP) :: cp
    real(DP) :: Massdisk, u_mdmtot
    real(DP) :: Luminosity

    real(DP) :: alpha(0:nar+1), Nu(0:nar+1), cs2(0:nar+1)
    real(DP) :: alpMRIeff(0:nar+1), alpGI(0:nar+1), alpdead(0:nar+1)
    real(DP) :: Ome(0:nar+1)
    real(DP) :: Md0, Md1, Md2, Mdaver
    real(DP) :: Mdr(0:nar+1)

    !sq=sqrt(), oo=1/()
    real(DP) :: sqR(0:nar+1), ooR(0:nar+1)
    real(DP) :: dR(1:nar+1), Rh(1:nar+1)
    real(DP) :: sqRh(1:nar+1), oodR(1:nar+1), oodRh(1:nar+1)

    integer :: it
  end type

  type(disk) :: md

  real(DP) :: Mstar, Lstar, Lacc

  !temporary variables
  real(DP) :: Sigm_temp, Sigg_temp, Q_temp, Te_temp
  real(DP) :: alpMRI_temp, alpdead_temp, alpGI_temp

contains

  subroutine init_star(Ms)
    implicit none
    real(DP), intent(in) :: Ms
    !set the initial mass of the Star
    !Mstar=Mstar_input*Msun
    Mstar=Ms
  end subroutine

  function luminosity_star(Starmass)
    implicit none
    real(DP) :: luminosity_star
    real(DP), intent(in) :: Starmass
    real(DP) :: lg10lum
    lg10lum=0.20_DP + 1.74_DP*log10(Starmass/Msun)
    luminosity_star=Lsun*10.0_DP**lg10lum
  end function

  subroutine init_disk(t)
    implicit none
    integer :: i
    real(DP) :: dx
    real(DP), intent(in) :: t

    !set the inital conditions

    md%t=t

    md%Mu=2.3_DP
    md%cp=2.7_DP*Rgas/md%Mu

    !inner and outer boundary conditions
    md%Rin=Rin_input  !cm 5Rsun
    md%Rout=Rout_input !cm 40AU

    dx=(sqrt(md%Rout)-sqrt(md%Rin))/dble(nar-1)
    !dx=(log(md%Rout)-log(md%Rin))/dble(nar-1)
    do i=0, nar+1
      md%R(i)=(dx*dble(i-1)+sqrt(md%Rin))**2
      !md%R(i)=exp(dble(i-1)*dx+log(md%Rin))
      md%ooR(i)=1.0_DP/md%R(i)
      md%sqR(i)=sqrt(md%R(i))
    end do

    do i=1, nar+1
      md%dR(i)=md%R(i)-md%R(i-1)
      md%oodR(i)=1.0_DP/md%dR(i)
      md%Rh(i)=(md%R(i)+md%R(i-1))/2.0_DP
      md%sqRh(i)=sqrt(md%Rh(i))
    end do

    do i=1, nar
      md%oodRh(i)=1.0_DP/(md%Rh(i+1)-md%Rh(i))
    end do

    !initial power low of Sigma and Tc
    do i=1, nar+1
      !md%Sig(i)=1200.0_DP*(md%R(i)/AU)**(-3.0_DP/2.0_DP)
      md%Sig(i)=0.0_DP
      md%Tc(i)=280.0_DP*(md%R(i)/AU)**(-1.0_DP/2.0_DP)
    end do

    !alpha_MRI
    do i=1, nar+1
      md%alpha(i)=alpha_MRI
    end do
  end subroutine

  subroutine boundary_cond()
    implicit none
    md%Sig(0)=0.0_DP
    !md%Sig(nar+1)=0.0_DP
    md%Tc(0)=20.0_DP
    md%Tc(nar+1)=1.0E-5_DP
  end subroutine

  subroutine explicit_evo(tdisk, dt)
    !explict solver of the evolution equation of disk
    implicit none
    type(disk), intent(inout) :: tdisk
    real(DP), intent(in) :: dt
    integer :: i

    real(DP) :: dSdt(nar), nsr(0:nar+1)
    real(DP) :: rdnsrdr(nar+1)

    do i=0, nar+1
      nsr(i)=tdisk%Nu(i)*tdisk%Sig(i)*tdisk%sqR(i)
    end do

    do i=1, nar+1
      rdnsrdr(i)=tdisk%sqRh(i)*(nsr(i)-nsr(i-1))*tdisk%oodR(i)
      tdisk%Mdr(i)=-6.0_DP*Pi*rdnsrdr(i)
    end do

    do i=1, nar
      dSdt(i)=3.0_DP*tdisk%ooR(i)*tdisk%oodRh(i)*(rdnsrdr(i+1)-rdnsrdr(i))
    end do


    do i=1, nar
      tdisk%Sig(i)=tdisk%Sig(i)+dSdt(i)*dt
    end do

    !tdisk%Sig(nar+1)=tdisk%Sig(nar)*tdisk%Nu(nar)*tdisk%sqR(nar)/tdisk%Nu(nar+1)/tdisk%sqR(nar+1)
    tdisk%Sig(nar+1)=tdisk%Sig(nar+1)-3.0_DP*tdisk%ooR(nar+1)*tdisk%oodRh(nar)*rdnsrdr(nar+1)*dt

    do i=1, nar
      tdisk%Sig(i)=max(tdisk%Sig(i), 1.0_DP)
    end do

  end subroutine

  subroutine implicit_evo(tdisk, dt)
    !implicit solver of the evolution equation of disk
    implicit none
    type(disk), intent(inout) :: tdisk
    real(DP), intent(in) :: dt
    integer :: i
    integer, save :: ifirst=0

    real(DP), save :: A(1:nar), B(1:nar), C(1:nar), BC(1:nar), D(1:nar)
    real(DP) :: up(1:nar-1), diag(nar), down(2:nar), u(nar)

    if (ifirst/=666) then
      do i=1, nar
        A(i)=-3.0_DP*tdisk%ooR(i)*tdisk%oodRh(i)*tdisk%sqRh(i+1)*tdisk%sqR(i+1)*tdisk%oodR(i+1)
        B(i)=3.0_DP*tdisk%ooR(i)*tdisk%oodRh(i)*tdisk%sqRh(i+1)*tdisk%sqR(i)*tdisk%oodR(i+1)
        C(i)=3.0_DP*tdisk%ooR(i)*tdisk%oodRh(i)*tdisk%sqRh(i)*tdisk%sqR(i)*tdisk%oodR(i)
        BC(i)=B(i)+C(i)
        D(i)=-3.0_DP*tdisk%ooR(i)*tdisk%oodRh(i)*tdisk%sqRh(i)*tdisk%sqR(i-1)*tdisk%oodR(i)
      end do
      ifirst=666
    end if

    do i=1, nar-1
      up(i)=dt*A(i)*tdisk%Nu(i+1)
    end do
    do i=2, nar
      down(i)=dt*D(i)*tdisk%Nu(i-1)
    end do
    do i=1, nar
      diag(i)=1.0_DP+dt*BC(i)*tdisk%Nu(i)
    end do

    call tridag_ser(down(2:nar),diag(1:nar), &
      up(1:nar-1),tdisk%Sig(1:nar),u(1:nar))

    !set boundary condition
    !no need

    do i=1, nar
      tdisk%Sig(i)=u(i)
    end do

    tdisk%Sig(nar+1)=tdisk%Sig(nar)*tdisk%Nu(nar)*tdisk%sqR(nar)/tdisk%Nu(nar+1)/tdisk%sqR(nar+1)

  end subroutine


  subroutine CN_evo(tdisk, dt)
    !CN solver of the evolution equation of disk
    !!!!NOT FINISHED!!!!
    !!!!DO NOT USE THIS!!!!
    implicit none
    type(disk), intent(inout) :: tdisk
    real(DP), intent(in) :: dt
    integer :: i
    integer, save :: ifirst=0
    real(DP) :: hdt

    real(DP), save :: A(1:nar), B(1:nar), C(1:nar), BC(1:nar), D(1:nar)
    real(DP) :: up(1:nar-1), diag(nar), down(2:nar), u(nar), rhs(nar)

    if (ifirst/=666) then
      do i=1, nar
        A(i)=-3.0_DP*tdisk%ooR(i)*tdisk%oodRh(i)*tdisk%sqRh(i+1)*tdisk%sqR(i+1)*tdisk%oodR(i+1)
        B(i)=3.0_DP*tdisk%ooR(i)*tdisk%oodRh(i)*tdisk%sqRh(i+1)*tdisk%sqR(i)*tdisk%oodR(i+1)
        C(i)=3.0_DP*tdisk%ooR(i)*tdisk%oodRh(i)*tdisk%sqRh(i)*tdisk%sqR(i)*tdisk%oodR(i)
        BC(i)=B(i)+C(i)
        D(i)=-3.0_DP*tdisk%ooR(i)*tdisk%oodRh(i)*tdisk%sqRh(i)*tdisk%sqR(i-1)*tdisk%oodR(i)
      end do
      ifirst=666
    end if

    hdt=dt/2.0_DP

    do i=1, nar-1
      up(i)=hdt*A(i)*tdisk%Nu(i+1)
    end do
    do i=2, nar
      down(i)=hdt*D(i)*tdisk%Nu(i-1)
    end do
    do i=1, nar
      diag(i)=1.0_DP+hdt*BC(i)*tdisk%Nu(i)
    end do

    rhs(1)=(-down(1+1)-diag(1)+2.0_DP)*tdisk%Sig(1)
    do i=2, nar-1
      rhs(i)=(-down(i+1)-diag(i)+2.0_DP-up(i-1))*tdisk%Sig(i)
    end do
    rhs(nar)=(-diag(nar)+2.0_DP-up(nar-1))*tdisk%Sig(nar)

    !do i=1, nar
    !  write(*,'(4G12.2)') down(i), diag(i), up(i), rhs(i)
    !end do
    !read(*,*)

    call tridag_ser(down(2:nar),diag(1:nar), &
      up(1:nar-1),rhs(1:nar),u(1:nar))

    !set boundary condition
    !no need

    do i=1, nar
      tdisk%Sig(i)=u(i)
    end do

  end subroutine

  subroutine add_Sigdot(tdisk, dt)
    !adding the mass infall to the disk
    implicit none
    type(disk), intent(inout) :: tdisk
    real(DP), intent(in) :: dt
    integer :: i

    do i=1, nar
      !tdisk%Sig(i)=tdisk%Sig(i)+Sigdotfun(tdisk%R(i))*dt
      tdisk%Sig(i)=tdisk%Sig(i)+Srtj_Core(tdisk%R(i), tdisk%t, tdisk%it)*dt
    end do
  end subroutine

  subroutine dt_fromTc(tdisk, dt, dTcdt)
    !calculating Qp and Qm
    !calculating dtmin
    implicit none
    type(disk), intent(inout) :: tdisk
    real(DP), intent(out) :: dt
    real(DP), intent(out) :: dTcdt(nar)
    integer :: i
    real(DP) :: Qp(nar), Qm(nar), dtft

    do i=1, nar
      !Qp(i)=Qp_fun(tdisk%Nu(i), tdisk%Sig(i), tdisk%Ome(i))
      Qp(i)=Qp_fun_allheating(tdisk, i)
      Qm(i)=Qm_fun(tdisk%Ome(i), tdisk%Tc(i), tdisk%Sig(i), tdisk%Mu)
      tdisk%Te(i)=Te_temp
      !Qm(i)=Qm_fun_Te(tdisk%Te(i))
      dTcdt(i)=2.0_DP*(Qp(i)-Qm(i))/tdisk%cp/tdisk%Sig(i)
    end do

    dt=1.0E3_DP*syear
    do i=1, nar
      dtft=0.1_DP*abs(tdisk%Tc(i)/dTcdt(i))
      if (dt>dtft) then
        dt=dtft
      end if
    end do

  end subroutine

  subroutine evo_Tc(tdisk, dt, dTcdt)
    !adding dT to Tc
    implicit none
    type(disk), intent(inout) :: tdisk
    real(DP), intent(in) :: dt, dTcdt(nar)
    integer :: i
    do i=1, nar
      tdisk%Tc(i)=tdisk%Tc(i)+dTcdt(i)*dt
      tdisk%Tc(i)=max(tdisk%Tc(i), 20.0_DP)
    end do
    tdisk%Tc(nar+1)=tdisk%Tc(nar)
  end subroutine

  subroutine calc_disk_acc(tdisk)
    !get disk accretion rate
    implicit none
    type(disk), intent(inout) :: tdisk
    real(DP) :: nsr1, nsr2, nsr3
    integer :: i

    nsr1=tdisk%Nu(1)*tdisk%Sig(1)*tdisk%sqR(1)
    nsr2=tdisk%Nu(2)*tdisk%Sig(2)*tdisk%sqR(2)
    nsr3=tdisk%Nu(3)*tdisk%Sig(3)*tdisk%sqR(3)
    tdisk%Md0=-6.0_DP*Pi*tdisk%sqRh(1)*(nsr1-0.0_DP)/tdisk%dR(1)
    tdisk%Md1=-6.0_DP*Pi*tdisk%sqRh(2)*(nsr2-nsr1)/tdisk%dR(2)
    tdisk%Md2=-6.0_DP*Pi*tdisk%sqRh(3)*(nsr3-nsr2)/tdisk%dR(3)

    tdisk%Mdaver=0.0_DP
    do i=iMds, iMde
      tdisk%Mdaver=tdisk%Mdaver+tdisk%Mdr(i)
    end do
    tdisk%Mdaver=tdisk%Mdaver/dble(iMde-iMds+1)
  end subroutine

  subroutine calc_disk_luminosity(tdisk)
    !get disk luminosity from Teff
    implicit none
    type(disk), intent(inout) :: tdisk
    integer :: i

    tdisk%Luminosity=0.0_DP
    do i=1, nar
      tdisk%Luminosity=tdisk%Luminosity+2.0_DP*Pi*tdisk%R(i)*tdisk%dR(i)*Stfb*tdisk%Te(i)**4
    end do
  end subroutine

  function acc_luminosity(mar, Starmass)
    implicit none
    !calculate disk accretion luminosity
    real(DP) :: acc_luminosity
    !mass accretion rate
    real(DP), intent(in) :: mar, Starmass

    !Bae et al 2014
    acc_luminosity=0.5_DP * Grav*(-mar)*Starmass/Rsun
    !Hartmann et al. (2016) uses 0.8 instead of 0.5
    !Here, we use 0.5 as in Bae et al 2014
  end function

  subroutine evolve(thedisk, Ms, tend)
    !evolution of the disk
    implicit none
    real(DP), intent(in) :: tend
    type(disk), intent(inout) :: thedisk
    real(DP), intent(inout) :: Ms
    real(DP) :: dt
    real(DP) :: dTcdt(nar)
    integer :: i

    !main loop
    do while (thedisk%t<tend)

      thedisk%Massdisk=mdisk()
      thedisk%u_mdmtot=thedisk%Massdisk/(thedisk%Massdisk+Ms)

      !calculate Nu_tot
      do i=0, nar+1
        thedisk%Ome(i)=Omega_fun(thedisk%R(i), Ms)
        thedisk%cs2(i)=cs2_fun(thedisk%Tc(i), thedisk%Mu)

       if (thedisk%Sig(i)==0.0_DP) then
        thedisk%alpMRIeff(i)=0.0_DP
        thedisk%QToom(i)=1.0E6_DP
        thedisk%alpGI(i)=0.0_DP
        thedisk%alpdead(i)=0.0_DP
        thedisk%alpha(i)=0.0_DP
        thedisk%Nu(i)=0.0_DP
       else

        !calc alp_MRI_eff, input(Tc, Sig), output(alpMRIeff, Sigm, Sigg=Sig-Sigm)
        call calc_alpMRIeff(thedisk%Tc(i), thedisk%Sig(i), thedisk%alpMRIeff(i), thedisk%Sigm(i), thedisk%Sigg(i))

        !calc QToom, then alpGI, using u_mdmtot if choose Kratter
        thedisk%QToom(i)=QToom(thedisk%cs2(i), thedisk%Ome(i), thedisk%Sig(i))
        thedisk%alpGI(i)=alpha_GI_Kratter(thedisk%QToom(i), thedisk%u_mdmtot)
        !thedisk%alpGI(i)=alpha_GI_Arimitage(thedisk%QToom(i))

        !thedisk%alpdead(i)=alpha_dead_const()
        thedisk%alpdead(i)=alpha_dead_Bae(thedisk%Sigm(i), thedisk%Sigg(i))

        !total alpha, Nu
        thedisk%alpha(i)=thedisk%alpMRIeff(i)+thedisk%alpGI(i)+thedisk%alpdead(i)
        thedisk%Nu(i)=nu_fun(thedisk%alpha(i), thedisk%cs2(i), thedisk%Ome(i))
       end if

        !call calc_nuet(thedisk%R(i), thedisk%Tc(i), thedisk%Sig(i), thedisk%alpha(i), thedisk%Mu, Ms, thedisk%Nu(i), thedisk%Te(i))
        !thedisk%Sigm(i)=Sigm_temp
        !thedisk%Sigg(i)=Sigg_temp
        !thedisk%QToom(i)=Q_temp
        !thedisk%alpMRIeff(i)=alpMRI_temp
        !thedisk%alpGI(i)=alpGI_temp
        !thedisk%alpdead(i)=alpdead_temp
        !thedisk%Nu(i)=nu_et(thedisk%R(i), thedisk%Tc(i), thedisk%Sig(i), thedisk%alpha(i), thedisk%Mu, Ms)
      end do

      call dt_fromTc(thedisk, dt, dTcdt)

      !set dt
      !dt=abs(thedisk%dR(i)/sqrt(thedisk%cs2(i)))
      !dt from evoTc
      do i=1, nar
        dt=min(dt, 0.02_DP*abs(thedisk%dR(i)**2/thedisk%Nu(i)))
      end do
      if (thedisk%t+dt>tend) then
        dt=tend-thedisk%t
      end if

      !viscous evolution
      call explicit_evo(thedisk, dt)
      !call implicit_evo(thedisk, dt)
      !call CN_evo(thedisk, dt)

      !adding mass infall
      call add_Sigdot(thedisk, dt)

      !Temperature evoution
      call evo_Tc(thedisk, dt, dTcdt)

      !output results of disk properties: Sigma, Tc, etc.
      call outfile(thedisk%t)

      !calculate disk accretion rate
      call calc_disk_acc(thedisk)
      !calculate disk luminosity
      call calc_disk_luminosity(thedisk)
      !calculate star luminosity
      Lstar=luminosity_star(Mstar)
      !calculate accretion luminosity
      Lacc=acc_luminosity(thedisk%Mdaver, Mstar)

      !add Mdot0 into Mstar
      Mstar=Mstar-md%Md0*dt
      !add Mass infall onto the Star from the core
      Mstar=Mstar+Mdot2Starfcore(thedisk%t)*dt

      call wriMtfile()

      thedisk%t=thedisk%t+dt
      thedisk%it=thedisk%it+1

      !test if Sigma is not a number
      !if nan, stop everything
      do i=1, nar
        if (isnan(thedisk%Sig(i))) then
          write(*,*) 'nan', i
          stop
        end if
      end do

    end do

    call outfile(thedisk%t)

  end subroutine

  subroutine outfile(t)
    !output disk properties
    implicit none
    real(DP), intent(in) :: t
    real(DP), save :: tnext
    character(15) :: fname, OutStr
    character(1) :: OldChar, NewChar
    integer :: i, j
    integer, save  :: ifirst=0
    !time interval for output
    real(DP), save :: tout

    tout=tout_input*syear

    if (t>tnext.OR.ifirst==0) then
      write(fname,'(A1, F10.1, A4)') 'a', t/syear, '.txt'

      OldChar=" "
      NewChar="0"
      OutStr=adjustl(trim(fname))
      j=index(OutStr,OldChar)
      do while(j>0)
        OutStr(j:j+LEN(OldChar)-1)=NewChar
        j=index(OutStr,OldChar)
      end do
      fname=OutStr
      write(*,*) fname
      open(88,file=fname)
      write(88, '(A4, F20.1, A4)') '# t=', t/syear, '  yr'
      write(88,'(A1, 7X, A1, A11, 7A20, 7X, A20)') &
        '#', 'i', 'R(AU)', 'Sigma', 'Tc', 'Te', &
        'Sigm', 'Sigg', 'Nu', 'Q', 'Mdot(Msun/yr)'
      do i=1, nar+1
        write(88,'(I10, 12G20.12)') i, &
          md%R(i)/AU, &
          md%Sig(i), &
          md%Tc(i), &
          md%Te(i), &
          md%Sigm(i), &
          md%Sigg(i), &
          md%Nu(i), &
          md%QToom(i), &
          md%Mdr(i)/Msun*syear, &
          md%alpMRIeff(i), &
          md%alpGI(i), &
          md%alpdead(i)
      end do
      close(88)
      tnext=tnext+tout
    end if

    if (ifirst/=666) then
      tnext=t+tout
      !write(*,*) 'tnext', tnext/syear
      ifirst=666
    end if

  end subroutine

  subroutine wriMtfile()
    !output the disk accretion rate
    implicit none
    integer, save :: ifirst=0
    real(DP), save :: Mdsave, Lsave, tnext

    if (ifirst/=666) then
      open(99,file='m.txt')
      write(99,'(A10, 5X, 3A20)') '#t(yr)', 'Mdot(Msun/yr)', 'Mdisk(Msun)', 'Mstar(Msun)'
      close(99)
      open(99,file='lumi.txt')
      write(99,'(A10, 5X, 4A20)') '#t(yr)', 'Lumi_disk(Lsun)', 'Lumi_acc(Lsun)', 'Lumi_star(Lsun)', 'Lumi_tot(Lsun)'
      close(99)
      Mdsave=md%Md2
      Lsave=md%Luminosity
      ifirst=666
      tnext=0.001_DP
    endif

    !if (abs((md%Md2-Mdsave)/Mdsave)>percent_Mdot_input) then
    !if (abs((md%Luminosity-Lsave)/Lsave)>percent_Mdot_input) then
    if (md%t > tnext) then
      open(99,file='m.txt',status='old',position='append')
      md%Massdisk=mdisk()
      write(99,'(4G20.12)') md%t/syear, -md%Mdaver/Msun*syear, md%Massdisk/Msun, Mstar/Msun
      close(99)
      open(99,file='lumi.txt',status='old',position='append')
      !md%Massdisk=mdisk()
      write(99,'(5G20.12)') md%t/syear, md%Luminosity/Lsun, Lacc/Lsun, Lstar/Lsun, (md%Luminosity+Lacc+Lstar)/Lsun
      close(99)
      Mdsave=md%Md2
      !Lsave=md%Luminosity
      tnext=md%t+syear*outMdot
    end if
  end subroutine

  function mdisk()
    implicit none
    real(DP) :: mdisk
    integer :: i
    mdisk=0.0_DP
    do i=1, nar
      mdisk=mdisk+2.0_DP*Pi*md%R(i)*md%dR(i)*md%Sig(i)
    end do
  end function

  function cs2_fun(T, mu)
    implicit none
    real(DP) :: cs2_fun
    real(DP), intent(in) :: T, mu
    cs2_fun=Rgas*T/mu
  end function

  function Omega_fun(R, Ms)
    implicit none
    real(DP) :: Omega_fun
    real(DP), intent(in) :: R, Ms
    Omega_fun=sqrt(Grav*Ms/R**3)
  end function

  function nu_fun(alp, cs2, Ome)
    implicit none
    real(DP) :: nu_fun
    real(DP), intent(in) :: alp, cs2, Ome
    nu_fun=alp*cs2/Ome
  end function

  function QToom(cs2, Ome, Sig)
    implicit none
    real(DP) :: QToom
    real(DP), intent(in) :: cs2, Ome, Sig
    QToom=sqrt(cs2)*Ome/Pi/Grav/Sig
  end function

  function alpha_GI_Zhu(QR)
    implicit none
    real(DP) :: alpha_GI_Zhu
    real(DP), intent(in) :: QR
    alpha_GI_Zhu=exp(-QR**4)
  end function

  function alpha_GI_Kratter(QR, u)
    !Kratter 2008 681:375 Eq.15-17 Q=max(Q,1) u=Md/(Md+Mstar) Eq.10
    implicit none
    real(DP) :: alpha_GI_Kratter
    real(DP), intent(in) :: QR, u
    real(DP) :: alp_short, alp_long, QR2
    QR2=max(1.0_DP, QR)
    alp_short=max(0.14_DP*(1.3_DP**2/QR2**2-1.0_DP)*(1.0_DP-u)**1.15_DP, 0.0_DP)
    alp_long=max((1.4E-3_DP*(2.0_DP-QR2))/(u**1.25_DP*sqrt(QR2)), 0.0_DP)
    alpha_GI_Kratter=sqrt(alp_short**2+alp_long**2)
    if (alpha_GI_Kratter>0.05_DP) then
      write(*,*) alpha_GI_Kratter
      alpha_GI_Kratter=0.05_DP
    end if
  end function

  function alpha_GI_Arimitage(Q)
    implicit none
    real(DP) :: alpha_GI_Arimitage
    real(DP), intent(in) :: Q
    if (Q<Qcirt_input) then
      alpha_GI_Arimitage=alpha_GI0*((Qcirt_input/Q)**2-1.0_DP)
    else
      alpha_GI_Arimitage=0.0_DP
    end if
    if (alpha_GI_Arimitage>0.05_DP) then
      write(*,*) alpha_GI_Arimitage
      alpha_GI_Arimitage=0.05_DP
    end if
  end function


  subroutine calc_alpMRIeff(T, Sig, alpMRIeff, Sig_a, Sig_d)
    implicit none
    real(DP), intent(out) :: alpMRIeff, Sig_a, Sig_d
    real(DP), intent(in) :: T, Sig

    if (T>Tcrit_input.OR.Sig<Sigcrit_input) then
      Sig_a=Sig
      Sig_d=0.0_DP
    else
      Sig_a=Sigcrit_input
      Sig_d=Sig-Sig_a
    end if

    alpMRIeff=Sig_a/Sig*alpha_MRI
  end subroutine

  function alpha_dead_const()
    implicit none
    real(DP) :: alpha_dead_const
    alpha_dead_const=alpha_Dead
  end function

  function alpha_dead_Bae(Sig_a, Sig_d)
    implicit none
    real(DP) :: alpha_dead_Bae
    real(DP), intent(in) :: Sig_a, Sig_d
    real(DP) :: frd
    frd=0.1_DP
    alpha_dead_Bae=min(1.0E-4_DP, frd*alpha_MRI*Sig_a/Sig_d)
  end function

  subroutine calc_nuet(R, T, Sig, alp, Mu, Ms, nu_et, Te)
    !effective viscosity Nu, the total viscosity
    implicit none
    real(DP), intent(out) :: nu_et, Te
    real(DP), intent(in) :: R, T, Sig, alp, Mu, Ms
    real(DP) :: nu_a, Sig_a, alp_G, nu_G, nu_D, Q, cs2, Ome
    real(DP) :: Sigm, Sigg, Tm, cs2m
    real(DP), save :: Sig_crit
    real(DP), save :: T_crit
    real(DP), save :: Q_crit

    Sig_crit=Sigcrit_input
    T_crit=Tcrit_input
    Q_crit=Qcirt_input

    if (Sig==0.0_DP) then
      nu_et=0.0_DP
      return
    end if

    !Sigm=sigact_fun(alp, T, Mstar, R, Sig)
    !Sigg=Sig-Sigm
    !Sigm_temp=Sigm
    !Sigg_temp=Sigg

    !Ome=Omega_fun(R, Ms)
    !call calcTmTe(alp, Sigm, Sigg, T, Ome, Mu, Tm, Te)
    Te=0.0_DP !just for out

    if (T>T_crit.OR.Sig<Sig_crit) then
      Sig_a=Sig
    else
      Sig_a=Sig_crit
    end if
    Sigm_temp=Sig_a
    Sigg_temp=Sig-Sig_a
    cs2=Rgas*T/Mu
    Ome=Omega_fun(R, Ms)
    nu_a=Sig_a/Sig*alp*cs2/Ome
    alpMRI_temp=Sig_a/Sig*alp
    !cs2m=Rgas*Tm/Mu
    !nu_a=Sigm/Sig*alp*cs2m/Ome


    Q=sqrt(cs2)*Ome/Pi/Grav/Sig
    Q_temp=Q
    if (Q<Q_crit) then
      alp_G=alpha_GI0*((Q_crit/Q)**2-1.0_DP)
    else
      alp_G=0.0_DP
    end if
    if (alp_G>0.05_DP) then
      write(*,*) alp_G
      alp_G=0.05_DP
    end if
    alpGI_temp=alp_G
    !nu_G=alp_G*cs2/Ome*Sigg/Sig
    nu_G=alp_G*cs2/Ome

    nu_D=alpha_Dead*cs2/Ome
    alpdead_temp=alpha_Dead

    nu_et=nu_a+nu_G+nu_D
  end subroutine


  function Sigdotfun(R)
    implicit none
    real(DP) :: Sigdotfun
    real(DP), intent(in) :: R
    real(DP) :: wid, Minfall
    Minfall=Mdotinfall_input*Msun/syear
    wid=dRadd_input*AU
    if (abs(R-Radd_input*AU)<5.0_DP*wid) then
      Sigdotfun=Minfall/wid/sqrt(2.0_DP*Pi)*exp(-(R-Radd_input*AU)**2/wid**2/2.0_DP)
      Sigdotfun=Sigdotfun/2.0_DP/Pi/Radd_input/AU
    else
      Sigdotfun=0.0_DP
    endif
  end function

  function opacityRossS2003(Temp, rho)
    implicit none
    real(8) :: opacityRossS2003
    real(8), intent(in) :: Temp, rho
    real(8), save :: opacRoss(10000,111)
    integer :: iT, ir
    integer, save :: ifirst=0

    if (ifirst/=123456) then
      open(444,file='kR.bin',form='unformatted',status='old')
      do iT=1, 10000
        read(444) opacRoss(iT,:)
      end do
      close(444)
      !write(*,*)'Rosseland Opacity Readed'
      ifirst=123456
    end if

    iT=nint(Temp)
    ir=nint((log10(rho)+18.0_8)/0.1_8)+1

    if (iT<5) then
      iT=5
    else if (iT>=9999) then
      iT=9999
    end if

    if (ir<1) then
      ir=1
    else if (ir>110) then
      ir=110
    end if

    opacityRossS2003=opacRoss(iT, ir)
  end function

  function opacityPlanckS2003(Temp, rho)
    implicit none
    real(8) :: opacityPlanckS2003
    real(8), intent(in) :: Temp, rho
    real(8), save :: opacPlanck(10000,111)
    integer :: iT, ir
    integer, save :: ifirst=0

    if (ifirst/=123456) then
      open(445,file='kP.bin',form='unformatted',status='old')
      do iT=1, 10000
        read(445) opacPlanck(iT,:)
      end do
      close(445)
      !write(*,*)'Planck Opacity Readed'
      ifirst=123456
    end if

    iT=nint(Temp)
    ir=nint((log10(rho)+18.0_8)/0.1_8)+1

    if (iT<5) then
      iT=5
    else if (iT>=9999) then
      iT=9999
    end if

    if (ir<1) then
      ir=1
    else if (ir>110) then
      ir=110
    end if

    opacityPlanckS2003=opacPlanck(iT, ir)
  end function

  function kappa(T, rho)
    !opacity
    implicit none
    real(DP) :: kappa
    real(DP), intent(in) :: T, rho
    real(DP), save :: Timax(12)
    real(DP), save :: k0(12), a(12), b(12)
    integer, save :: ifirst=0
    integer :: i

    !Eq 10 of Armitage et al 2001
    !kappa=0.02_DP*T**0.8_DP
    !return

    !opacity table from Bell & Lin
    !see Table 1 of Armitage et al 2001
    if (ifirst/=666) then
      k0(1)=1.00E-04_DP
      k0(2)=3.00E+00_DP
      k0(3)=1.00E-02_DP
      k0(4)=5.00E+04_DP
      k0(5)=1.00E-01_DP
      k0(6)=2.00E+15_DP
      k0(7)=2.00E-02_DP
      k0(8)=2.00E+81_DP
      k0(9)=1.00E-08_DP
      k0(10)=1.00E-36_DP
      k0(11)=1.50E+20_DP
      k0(12)=3.48E-01_DP
      a(1)=0.00_DP
      a(2)=0.00_DP
      a(3)=0.00_DP
      a(4)=0.00_DP
      a(5)=0.00_DP
      a(6)=0.00_DP
      a(7)=0.00_DP
      a(8)=1.00_DP
      a(9)=0.6666666666666666667_DP
      a(10)=0.333333333333333333_DP
      a(11)=1.00_DP
      a(12)=0.00_DP
      b(1)=2.10_DP
      b(2)=-0.01_DP
      b(3)=1.10_DP
      b(4)=-1.50_DP
      b(5)=0.70_DP
      b(6)=-5.20_DP
      b(7)=0.80_DP
      b(8)=-24.00_DP
      b(9)=3.00_DP
      b(10)=10.00_DP
      b(11)=-2.50_DP
      b(12)=0.00_DP

      do i=1, 6
        Timax(i)=(k0(i)/k0(i+1))**(1.0_DP/(b(i+1)-b(i)))
      end do

      ifirst=666
    end if

    do i=7,12
      Timax(i)=(k0(i)/k0(i+1)*rho**(a(i)-a(i+1)))**(1.0_DP/(b(i+1)-b(i)))
    end do

    do i=1, 12
      if (T<=Timax(i)) then
        kappa=k0(i)*rho**a(i)*T**b(i)
        return
      end if
    end do
    kappa=k0(12)

  end function

  function Qp_fun(nu, Sig, Ome)
    !viscous heating
    implicit none
    real(DP) :: Qp_fun
    real(DP), intent(in) :: nu, Sig, Ome
    Qp_fun=9.0_DP/8.0_DP*nu*Sig*Ome**2
  end function

  function Qp_fun_allheating(td, i)
    implicit none
    type(disk), intent(inout) :: td
    real(DP) :: Qp_fun_allheating
    !viscous, infall, exteranl heating
    real(DP) :: Qp_vis, Qp_infall, Qp_ext
    !star, accretion, envelope heating
    real(DP) :: sT4_star, sT4_acc, sT4_env
    real(DP) :: Rc, RRc
    real(DP) :: fstar, facc
    real(DP) :: tau, H, kap, rho
    integer, intent(in) :: i
    !viscous heating
    Qp_vis=9.0_DP/8.0_DP*td%Nu(i)*td%Sig(i)*td%Ome(i)**2

    !infall heating Eq 19 of Bae et al 2013
    Rc=Rd_Core(td%t)
    RRC=td%R(i)/Rc
    if (RRc<0.2_DP) then
      Qp_infall=0.0_DP
    else
      Qp_infall=Grav*Mstar*Mdot_Core()/4.0_DP/Pi/Rc**3 *(3.0_DP-2.0_DP*sqrt(RRc))/RRc**2
    end if

    fstar=0.1_DP
    sT4_star=fstar*Lstar/4.0_DP/Pi/td%R(i)**2

    facc=0.1_DP
    sT4_acc=facc*Lacc/4.0_DP/Pi/td%R(i)**2

    sT4_env=Stfb*Temp_Core**4


    H=sqrt(Rgas*td%Tc(i)/td%Mu)/td%Ome(i)
    rho=td%Sig(i)/2.0_DP/H
    !kap=kappa(Tc, rho)
    kap=opacityRossS2003(td%Tc(i), rho) !Use Semenov et al 2003 opacity
    tau=kap*td%Sig(i)

    !Eq 12 of Bae et al 2013
    Qp_ext=16.0_DP/3.0_DP*(sT4_star+sT4_acc+sT4_env)*tau/(1.0_DP+tau**2)

    !Eq 22 of Bae et al 2013
    Qp_fun_allheating=Qp_vis+Qp_infall+Qp_ext

  end function

  function Qm_fun(Ome, Tc, Sig, Mu)
    !Qm of Armitage et al 2001
    implicit none
    real(DP) :: Qm_fun
    real(DP), intent(in) :: Ome, Tc, Sig, Mu
    real(DP) :: H, rho, Te4, kap, tau, Q

    H=sqrt(Rgas*Tc/Mu)/Ome
    rho=Sig/2.0_DP/H
    Q=sqrt(Rgas*Tc/Mu)*Ome/Pi/Grav/Sig
    !kap=kappa(Tc, rho)
    kap=opacityRossS2003(Tc, rho) !Use Semenov et al 2003 opacity
    if (Q<Qcirt_input) then
      tau=Sig/2.0_DP*kap
    else
      if (Tc>Tcrit_input .OR. Sig<Sigcrit_input) then
        tau=Sig/2.0_DP*kap
      else
        tau=Sigcrit_input/2.0_DP*kap
      endif
    end if

    !tau=Sig/2.0_DP*kap

    Te4=8.0_DP/3.0_DP/tau*Tc**4
    Te4=min(Te4, Tc**4)
    Te_temp=Te4**0.25_DP
    Qm_fun=Stfb*Te4
  end function

  function Qm_fun_Te(Te)
    implicit none
    real(DP), intent(in) :: Te
    real(DP) :: Qm_fun_Te
    Qm_fun_Te=Stfb*Te**4
  end function


  function sigact_fun(alp, T, Mstar, R, Sig)
    !sigma active for Martin et al 2012
    !Sigma_active can not be greater than sigma
    implicit none
    real(DP), intent(in) :: alp, T, Mstar, R, Sig
    real(DP) :: sigact_fun

    sigact_fun=min(sigact_martin(alp, T, Mstar, R), Sig)
    !sigact_fun=Sig !for fully active
  end function

  function sigact_martin(alp, T, Mstar, R)
    ! Sigma active, maximun of
    ! Sigma active for cosmic ray ionization
    ! and Sigma active for thermal ionization
    ! Eq(27) of Martin et al 2012, MNRAS, 420, 3239
    implicit none
    real(DP), intent(in) :: alp, T, Mstar, R
    real(DP) :: Sigact_martin
    Sigact_martin=max(Sigact_cr(alp, T, Mstar, R), &
      Sigact_th(alp, T, Mstar, R))
  end function

  function Sigact_cr(alp, T, Mstar, R)
    ! Sigma active for cosmic ray ionization
    ! Eq(27) of Martin et al 2012, MNRAS, 420, 3239
    implicit none
    real(DP), intent(in) :: alp, T, Mstar, R
    real(DP) :: Sigact_cr
    real(DP) :: Remcrit=1.0E4_DP
    real(DP) :: zeta=1.0E-17_DP !s^-1

    Sigact_cr=1.36_DP*(alp/0.1_DP)*(T/100.0_DP)**2 &
      *(Remcrit/1.0E4_DP)**(-2.0_DP) &
      *(Mstar/Msun)**(-3.0_DP/2.0_DP) &
      *(zeta/1.0E-17_DP) &
      *(R/AU)**(9.0_DP/2.0_DP)

  end function


  function Sigact_th(alp, T, Mstar, R)
    ! Sigma active for thermal ionization
    ! Eq(28) of Martin et al 2012, MNRAS, 420, 3239
    implicit none
    real(DP), intent(in) :: alp, T, Mstar, R
    real(DP) :: Sigact_th
    real(DP) :: Remcrit=1.0E4_DP
    real(DP) :: Tion=1500.0_DP

    Sigact_th=0.2_DP*(alp/0.1_DP)*(T/Tion)**3 &
      *(Remcrit/1.0E4_DP)**(-2.0_DP) &
      *(R/0.03_DP/AU)**(+9.0_DP/2.0_DP) &
      *(Mstar/Msun)**(-3.0_DP/2.0_DP) &
      *exp(-25188.0_DP*(1.0_DP/T-1.0_DP/Tion)*2.0_DP)

  end function

  subroutine calcTmTe(alp, sigm, sigg, Tc, Ome, Mu, Tm, Te)
    implicit none
    real(DP), intent(in) :: alp, sigm, sigg, Tc, Ome, Mu
    real(DP), intent(out) :: Tm, Te
    real(DP) :: taum, taug, taum1, taum2
    real(DP) :: Tm1, Te1, Tm2, Te2

    !if (sigg<1.0E-90_DP) then
    if (0>1) then
      !sigg=0 no dead zone, fully active
      !T_m=T_c
      Tm=Tc
      !calc tau_m
      !Eq(6) of Martin et al 2012
      taum=tau_fun(Tm, sigm, Ome, Mu)
      !write(*,*) 'sm, sg, taum', sigm, sigg, taum
      !write(*,*) 'Tm, Tc', Tm, Tc

      if (taum<1.0_DP) then
        !thin
        Te=Tm
        !Tc=Tm=Te
      else
        !thick
        Te=Tm/taum**0.25_DP
      end if

    else
      !sigg>0, dead zone
      !assuming thick
      taug=tau_fun(Tc, sigg, Ome, Mu)
      call calcTmTethick(alp, taug, Tc, sigm, Ome, Mu, Tm1, Te1, taum1)

      !assuming thin
      !Tm=Te
      !Eq(12) of Martin et al 2012
      Tm2=(Tc**4/taug)**0.25_DP
      Te2=Tm2
      taum2=tau_fun(Tm2, sigm, Ome, Mu)

      if (taum1>=1.0_DP) then
        !1 thick
        if (taum2>1.0_DP) then
          !2 thick
          !write(*,*) 'thick'
          Tm=Tm1
          Te=Te1
        else
          !2 thin
          !both thick and thin
          !write(*,*) 'thick thin'
          !stop
          Tm=Tm1
          Te=Te1
        end if
      else
        !1 thin
        if (taum2>1.0_DP) then
          !2 thick
          !no solution
          !write(*,*) 'neither thick nor thin'
          !stop
          Tm=Tm1
          Te=Te1
        else
          !2 thin
          ! thin
          !write(*,*) 'thin'
          Tm=Tm2
          Te=Te2
        end if
      end if

    end if

  end subroutine

  subroutine calcTmTethick(alp, taug, Tc, sigm, Ome, Mu, Tm, Te, taum)
    !simple interation method for Tm
    !the get Te
    implicit none
    real(DP), intent(in) :: alp, taug, Tc, sigm, Ome, Mu
    real(DP), intent(out) :: Tm, Te, taum
    integer :: iter
    real(DP) :: Tmn, Tmo

    iter=1
    Tmo=Tc
    do
      Tmn=Tmfun(alp, taug, Tc, sigm, Ome, Tmo, Mu)
      iter=iter+1
      if (iter>30) then
        write(*,*) 'iter>30'
        stop
      end if
      if (abs(Tmn-Tmo)/Tmo>0.001_DP) then
        Tmo=Tmn
      else
        Tm=Tmn
        taum=tau_fun(Tm, sigm, Ome, Mu)
        Te=(Tm**4/taum)**0.25_DP
        exit
      end if
    end do

  end subroutine

  function Tmfun(alp, taug, Tc, sigm, Ome, Tm, Mu)
    !Eq(11) of Martin et al 2012
    implicit none
    real(DP), intent(in) :: alp, taug, Tc, sigm, Ome, Tm, Mu
    real(DP) :: Tmfun
    real(DP) :: num, cs2m, Tm4, taum
    cs2m=cs2_fun(Tm, Mu)
    num=alp*cs2m/Ome
    taum=tau_fun(Tm, sigm, Ome, Mu)
    Tm4=(Tc**4+9.0_DP/8.0_DP/Stfb*num*sigm*Ome**2*taug)/(1.0_DP+taug/taum)
    Tmfun=Tm4**0.25_DP
  end function

  function tau_fun(T, sig, Ome, Mu)
    !Eq(6,7) of Martin et al 2012
    implicit none
    real(DP), intent(in) :: T, sig, Ome, Mu
    real(DP) :: tau_fun
    real(DP) :: rho
    rho=sig/2.0_DP/sqrt(cs2_fun(T, Mu))*Ome
    tau_fun=3.0_DP/16.0_DP*kappa(T, rho)*sig
    !write(*,*) 'T, sig, kappa(T, rho), tau_fun'
    !write(*,*) T, sig, kappa(T, rho), tau_fun
  end function

end module



program evol
  use CloudCoreMod
  use StellarMod
  use evolve_mod
  implicit none
  real(DP) :: t
  integer :: it

  call get_parameters()

  call ReadCloudCorePara('cloudcorepara.txt')
  call ReadStellarPara('Stellarpara.txt')

!----------------------------------------------------------
  t=0.0_DP*syear
  do while (Rd_Core(t)<Rin_input*1.1_DP)
    t=t+1.0_DP*syear
  enddo
  write(*,*) 't start at', t/syear
  it=0
  Mass_Stellar=Mdot_Core()*t
  write(*,*) 'Mass_Stellar', Mass_Stellar/Msun
!----------------------------------------------------------

  md%it=0

  call init_star(Mass_Stellar)
  call init_disk(t)
  call boundary_cond()

  call evolve(md, Mstar, tend_input*syear)

end program
