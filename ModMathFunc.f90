!> @file ModMathFunc.f90
!> @brief Module for mathematical functions
!> @author Masaaki Kimura
module ModMathFunc
  use iso_fortran_env, only: real64, real128, error_unit
  use ModNuclConst, only: PI, VDIM, VX, VY, VZ, SDIM, SU, SD, I_IMAG
  implicit none
  private

  public :: MathFunc !> public class for math functions


  ! derived type for math functions which is a collection of math functions
  type MathFunc 

    contains

    !>!> gamma function  !>!>
    procedure, nopass, private :: dgamln_  ! real(real64) version
    procedure, nopass, private :: qgamln_  ! real(real128) version
    generic, public :: gamln => dgamln_, qgamln_

    !>!> associated legendre function !>!>
    procedure, nopass, private :: dplgndr_  ! real(real64) version
    procedure, nopass, private :: zplgndr_  ! complex(real64) version
    generic, public :: plgndr => dplgndr_, zplgndr_

    ! Shperical Harmonics from polar angles    
    procedure, nopass, public :: sph_harm
    ! Shperical Harmonics from real/complex valued 3D vector
    generic, public :: sph_harm_vector => dsph_harm_, zsph_harm_
    procedure, nopass, private :: dsph_harm_ ! real-valued 3D vector
    procedure, nopass, private :: zsph_harm_ ! complex-valued 3D vector

    ! Wigner's D function
    procedure, nopass :: wigner_d => dwigner_d_

    ! Rotation matrix
    procedure, nopass :: rot_cart  ! calculates rotation matrix in cartesian rep.
    procedure, nopass :: rot_spin  ! calculates rotation matrix in spinor rep.
    procedure, nopass :: euler_from_rot_cart ! calculates euler angles from rotation matrix

    ! Clebsch-Gordan coefficient
    procedure, nopass :: cg_coeff
    
    ! gauss-legendre quadrature
    procedure, nopass, public :: gauleg
  end type MathFunc

  
  contains

  !> @brief returns logarithm of gamma function in double precision
  !> @returns logarithm of gamma function in double precision
  !> @param xx argument for gamma function
  pure elemental real(real64) function dgamln_(xx)
    implicit none
    ! arguments
    real(real64), intent(in) :: xx
    ! constants
    real(real64), parameter :: cof(6) = &
    [ 76.18009172947146_real64, -86.50532032941677_real64, 24.01409824083091_real64, &
    & -1.231739572450155_real64, 0.1208650973866179e-2_real64, -0.5395239384953e-5_real64 ]
    ! internal variables
    integer :: j
    real(real64) :: x, y, ser, tmp
    
    ! body
    x = xx
    y = x
    tmp = x + 5.5_real64
    tmp = tmp - (x + 0.5_real64)*log(tmp)
    ser = 1.000000000190015_real64

    do j=1, 6
      y = y + 1.0_real64
      ser = ser + cof(j)/y
    end do
    dgamln_ = -tmp + log(2.5066282746310005_real64*ser/x)
  end function dgamln_
  
  
  !> @brief returns logarithm of gamma function in quadratic precision
  !> @returns logarithm of gamma function in quadratic precision
  !> @param x argument for gamma function
  pure elemental real(real128) function qgamln_(x)
    implicit none
    ! arguments 
    real(real128), intent(in) :: x
    ! constants
    real(real128), parameter :: pv= 1.357812200070394647397691360527351888e01_real128
    real(real128), parameter :: pr= 3.1782384299734898421289539143998119380e-06_real128
    real(real128), parameter :: p0= 3.1482070283349300354582623623908339457e-04_real128
    real(real128), parameter :: p1= 1.2793741608722984500693458490473661860e-02_real128
    real(real128), parameter :: p2= 2.7874830306029980874434569055259616606e-01_real128
    real(real128), parameter :: p3= 3.5748763958228570180758258557929027134e+00_real128
    real(real128), parameter :: p4= 2.7927280421563325015666935178375281217e+01_real128
    real(real128), parameter :: p5= 1.3321384650379738989446885832268784755e+02_real128
    real(real128), parameter :: p6= 3.7950405192465422312792634449147935784e+02_real128
    real(real128), parameter :: p7= 6.1562149993028259463346808196235292341e+02_real128
    real(real128), parameter :: p8= 5.2400400869100650701118261358974985117e+02_real128
    real(real128), parameter :: p9= 2.0418766202023711876168179075996496480e+02_real128
    real(real128), parameter :: p10= 2.8645619772729108683191342647193554200e+01_real128
    real(real128), parameter :: p11= 8.9507210141338984737305834751291040398e-01_real128
    real(real128), parameter :: p12= 1.8410863315761265630602733481713520754e-03_real128
    real(real128), parameter :: q1= 9.9999999999999999999999999982917706744e-01_real128
    real(real128), parameter :: q2= 2.0000000000000000000000014372532510977e+00_real128
    real(real128), parameter :: q3= 2.9999999999999999999987806987019196983e+00_real128
    real(real128), parameter :: q4= 4.0000000000000000003182289605630538899e+00_real128
    real(real128), parameter :: q5= 4.9999999999999999603276569479688692879e+00_real128
    real(real128), parameter :: q6= 6.0000000000000030034309156698097129604e+00_real128
    real(real128), parameter :: q7= 6.9999999999998375247462698288225315906e+00_real128
    real(real128), parameter :: q8= 8.0000000000071915518803021765161684809e+00_real128
    real(real128), parameter :: q9= 8.9999999997000681861822653951282600849e+00_real128
    real(real128), parameter :: q10= 1.0000000014205005237309132429530491661e+01_real128
    real(real128), parameter :: q11= 1.0999998953920119680361242478373085334e+01_real128
    real(real128), parameter :: q12= 1.2000238108934194337280539725944422661e+01_real128
    ! local variables
    real(real128) :: w,y 

    ! body
    w=x
    if(x.lt.0) w=1-x
    y=log(((((((((((((p12/(w+q12)+p11)/(w+q11)+p10)/(w+q10)+p9)/(w+q9)+p8)/(w+q8)+p7)/ &
      & (w+q7)+p6)/(w+q6)+p5)/(w+q5)+p4)/(w+q4)+p3)/(w+q3)+p2)/(w+q2)+p1)/(w+q1)+p0)/  &
      &  w+pr)+(w-0.5_real128)*log(w+pv)-w
    if(x.lt.0) y=log(pi/sin(pi*x))-y
    qgamln_=y
  end function qgamln_
  
  !> @brief real(real64) version of the Legendre polynominals\
  !> @note l and m are not doubled as in wigner's D function and clebsch-gordan coefficient
  !> @returns associated Legendre polynominals P_l^m(x) in real(real64) precision
  !> @param l rank of the polynomial
  !> @param m order of the polynomial (0 <= m <= l)
  !> @param x argument of the polynomial
  pure elemental real(real64) function dplgndr_(l,m,x)
    implicit none
    ! argument
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(real64), intent(in) :: x
    ! internal
    integer :: i,ll,am
    real(real64) :: fact,pll,pmm,pmmp1,somx2
    
    pll = 0 ! to suppress compiler warning
    am = abs(m)
    pmm = 1.0_real64
    if ( am > 0) then
      somx2 = sqrt(1.0d0 - x*x)
      fact = 1.0_real64
      do i=1,am
        pmm = -pmm*fact*somx2
        fact = fact + 2.0_real64
      end do
    end if
    if(l==am) then
      dplgndr_ = pmm
    else
      pmmp1 = x*(2*am+1)*pmm
      if(l == am+1) then
        dplgndr_ = pmmp1
      else
        do ll=am+2,l
          pll = (x*(2*ll-1)*pmmp1-(ll+am-1)*pmm)/(ll-am)
          pmm = pmmp1
          pmmp1 = pll
        end do
        dplgndr_ = pll
      end if
    end if
    
    if(m < 0) dplgndr_ = dplgndr_*(-1)**am
  end function dplgndr_

  !> @brief complex(real64) version of the associated Legendre polynomials\
  !> @note l and m are not doubled as in wigner's D function and clebsch-gordan coefficient
  !> @returns associated Legendre polynomials P_l^m(z) in complex(real64) precision
  !> @param l rank of the polynomial
  !> @param m order of the polynomial (0 <= m <= l)
  !> @param z argument of the polynomial in complex(real64) precision
  pure elemental complex(real64) function zplgndr_(l,m,z)
    implicit none
    ! arguments
    integer, intent(in) :: l,m
    complex(real64), intent(in) :: z
    ! Computes the associated Legendre polynomial Pml (z). 
    ! Here m and l are integers satisfying 0 ≤ m ≤ l, while z lies in the range −1 ≤ |z| ≤ 1.
    ! local variables
    integer :: i,ll
    real(real64) :: fact
    complex(real64) :: pmm,pmmp1,somx2

    ! body
    ! Compute Pmm which has a simle analytic form
    zplgndr_ = 0.0_real64 ! to suppress compiler warning
    pmm = 1.0_real64
    somx2 = sqrt((1.0_real64-z)*(1.0_real64+z))
    fact=1.0_real64
    do i=1,m
      pmm = -pmm*fact*somx2
      fact = fact+2.0_real64
    end do

    ! Compute Plm using recurrence formula (See numerical recipe 'plgndr')
    if (l == m) then ! l == m case
      zplgndr_ = pmm
      return
    end if

    pmmp1 = z*(2*m+1)*pmm
    if(l == m+1) then ! l == m+1 case
      zplgndr_ = pmmp1
      return
    end if

    do ll=m+2,l ! l > m+1 cases, we use recurrence formula
      zplgndr_=(z*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
      pmm=pmmp1
      pmmp1=zplgndr_
    end do
  end function zplgndr_


  !>
  !> @brief spherical harmonics from polar angles\
  !> @note l and m are not doubled as in wigner's D function and clebsch-gordan coefficient
  !> @returns spherical harmonics Y_l^m(theta,phi) in complex(real64) precision
  !> @param l orbital angular momentum
  !> @param m magnetic substitute
  !> @param theta polar angle
  !> @param phi azimuthal angle
  pure elemental complex(real64) function sph_harm(l,m,theta,phi)
    implicit none
    ! argument
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(real64), intent(in) :: theta
    real(real64), intent(in) :: phi
    ! internal
    integer :: i,ll,am
    real(real64) :: fact,pll,pmm,pmmp1,somx2,x,cnst
    complex(real64) :: e
    
    pll = 0.0_real64 ! to suppress compiler warning
    am = abs(m)
    x = cos(theta)
    e = exp(cmplx(0.0_real64,am*phi,kind=real64))
    cnst = sqrt((2*l+1)/(4*PI)*exp(dble(qgamln_(l-am+1.0_real128)-qgamln_(l+am+1.0_real128))))
    
    pmm = 1.0_real64
    if ( am > 0) then
      somx2 = sqrt(1.0_real64 - x*x)
      fact = 1.0_real64
      do i=1,am
        pmm = -pmm*fact*somx2
        fact = fact + 2.0_real64
      end do
    end if
    if(l==am) then
      sph_harm = pmm
    else
      pmmp1 = x*(2*am+1)*pmm
      if(l == am+1) then
        sph_harm = pmmp1
      else
        do ll=am+2,l
          pll = (x*(2*ll-1)*pmmp1-(ll+am-1)*pmm)/(ll-am)
          pmm = pmmp1
          pmmp1 = pll
        end do
        sph_harm = pll
      end if
    end if
    sph_harm = cnst*e*sph_harm
    
    if(m < 0) sph_harm = conjg(sph_harm)*(-1)**am
    
  end function sph_harm

  !>
  !> @brief calculate spherical harmonics from real-valued 3D vector\
  !> @note l and m are not doubled as in wigner's D function and clebsch-gordan coefficient
  !> @returns spherical harmonics Y_l^m(r) in complex(real64) precision
  !> @param l orbital angular momentum
  !> @param m magnetic substitute
  !> @param r real-valued 3D vector
  pure complex(real64) function dsph_harm_(l, m, r)
    implicit none
    ! arguments
    integer, intent(in) :: l ! orbital angular momentum
    integer, intent(in) :: m ! magnetic substitute
    real(real64), intent(in) :: r(VDIM) ! real-valued three dimensional vector
    ! internal variables
    real(real64) :: rr, costheta, cnst
    complex(real64) :: xpy, xmy

    ! body
    ! return Y00 for l=m=0 case for all vectors
    if(l==0) then
      dsph_harm_ = 1.0_real64/sqrt(4*PI)
      return
    end if

    ! compute the 'length' of real-valued 3D vector r
    rr = sqrt(dot_product(r(:),r(:)))
    ! if the length is zero, return zero for l/=0 cases
    if (.not. (rr > 0.0_real64)) then
      dsph_harm_ = 0.0_real64
      return
    end if

    costheta = r(VZ)/rr ! cos(theta) 
    if(.not. (r(VX)*r(VX) + r(VY)*r(VY) > 0.0_real64)) then
      xpy = 1.0_real64 ! Y function does not depend on phi if r*sin(theta) == 0.0d0
      xmy = 1.0_real64
    else
      xpy = (r(VX) + I_IMAG*r(VY))/sqrt(r(VX)*r(VX) + r(VY)*r(VY)) ! exp(i\phi)
      xmy = (r(VX) - I_IMAG*r(VY))/sqrt(r(VX)*r(VX) + r(VY)*r(VY)) ! exp(-i\phi)
    end if
    ! compute normalization factor
    cnst = sqrt((2*l+1)/(4*PI)*exp(dble(qgamln_(l-abs(m)+1.0_real128)-qgamln_(l+abs(m)+1.0_real128))))
    if (m >= 0) then
      dsph_harm_ = cnst*dplgndr_(l,abs(m),costheta)*xpy**abs(m)
    else
      dsph_harm_ = cnst*dplgndr_(l,abs(m),costheta)*(-xmy)**abs(m)
    end if

  end function dsph_harm_

  !>
  !> @brief spherical harmonics from complex-valued 3D vector\
  !> @note l and m are not doubled as in wigner's D function and clebsch-gordan coefficient
  !> @returns spherical harmonics Y_l^m(z) in complex(real64) precision
  !> @param l orbital angular momentum
  !> @param m magnetic substitute
  !> @param z complex-valued 3D vector
  pure complex(real64) function zsph_harm_(l, m, z)
    ! arguments
    integer, intent(in) :: l ! orbital angular momentum
    integer, intent(in) :: m ! magnetic substitute
    complex(real64), intent(in) :: z(VDIM) ! complex valued three dimensional vector
    ! internal variables
    real(real64) :: cnst
    complex(real64) :: r, costheta, xpy, xmy

    ! body
    ! return Y00 for l=m=0 case for all vectors
    if(l==0) then
      zsph_harm_ = 1.0_real64/sqrt(4*PI)
      return
    end if

    ! compute the 'length' of complex-valued 3D vector z
    r = sqrt(sum(z(:)*z(:))) ! this is sqrt of z.z (NOT z*.z)
    ! if the length is zero, return zero.
    if (.not. (abs(real(r,real64)) + abs(aimag(r)) > 0.0_real64)) then
      zsph_harm_ = 0.0_real64
      return
    end if

    costheta = z(VZ)/r ! this corresponds to cos(theta) for the real vector case
    if(.not. (abs(real(z(VX)*z(VX) + z(VY)*z(VY),real64)) + abs(aimag(z(VX)*z(VX) + z(VY)*z(VY))) > 0.0_real64)) then
      xpy = 1.0_real64 ! Y function does not depend on phi if r*sin(theta) == 0.0d0
      xmy = 1.0_real64
    else
      xpy = (z(VX) + I_IMAG*z(VY))/sqrt(z(VX)*z(VX) + z(VY)*z(VY)) ! this corresponds to exp(i\phi)
      xmy = (z(VX) - I_IMAG*z(VY))/sqrt(z(VX)*z(VX) + z(VY)*z(VY)) ! this corresponds to exp(-i\phi)
    end if
    ! compute normalization factor
    cnst = sqrt((2*l+1)/(4*PI)*exp(dble(qgamln_(l-abs(m)+1.0_real128)-qgamln_(l+abs(m)+1.0_real128))))
    if (m >= 0) then
      zsph_harm_ = cnst*zplgndr_(l,abs(m),costheta)*xpy**abs(m)
    else
      zsph_harm_ = cnst*zplgndr_(l,abs(m),costheta)*(-xmy)**abs(m)
    end if

  end function zsph_harm_

  
  
  !>
  !> @brief Wigner's D function\
  !> @note j, m1 and m2 are doubled (j=1 means angular momentum 1/2)
  !> @returns Wigner's D function D^j_{m1,m2}(alpha,beta,gamma) in complex(real64) precision
  !> @param j rank of the D function
  !> @param m1 first magnetic substitute
  !> @param m2 second magnetic substitute
  !> @param alp alpha of Euler angles
  !> @param bet beta of Euler angles
  !> @param gam gamma of Euler angles
  pure function dwigner_d_(j,m1,m2,alp,bet,gam)
    implicit none
    ! arguments
    integer, intent(in) :: j
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    real(real64), intent(in) :: alp(:)
    real(real64), intent(in) :: bet(:)
    real(real64), intent(in) :: gam(:)
    ! return
    complex(real64), allocatable :: dwigner_d_(:,:,:)
    ! internal variables
    integer :: k,mink,maxk,step,A,B,C,kapp,ia,ig,ib,asz,gsz,bsz
    real(real64) :: cb,sb,cnst
    real(real64) :: fct(max((j+m1)/2,(j-m1)/2,(j+m2)/2,(j-m2)/2)+1)
    real(real64), allocatable :: coff(:), db(:)
    complex(real64), allocatable :: ea(:), eg(:)
    
    ! body
    ! number of alpha, gamma and beta grids
    asz = size(alp)
    bsz = size(bet)
    gsz = size(gam)
    allocate(dwigner_d_(asz,bsz,gsz),db(bsz))
    
    ! calculate ln(factorial), fct(n) = ln( (n-1)! )
    fct(:) = 0.0_real128
    do k=3, size(fct(:))
      fct(k) = dble(qgamln_(k+0.0_real128))
    end do

    ! calculate exp(-im1*alp) and exp(-im2*gam)
    ea = exp(-0.50_real64*m1*alp(:)*I_IMAG)
    eg = exp(-0.50_real64*m2*gam(:)*I_IMAG)

    ! calculate coefficients of small d function
    A = (j-m1)/2
    B = (j+m2)/2
    C = (m1-m2)/2
    cnst = (fct((j+m1)/2+1) + fct((j-m2)/2+1) + fct(A+1) + fct(B+1))/2
    mink = max(0,-C)
    maxk = min(A,B)
    allocate(coff(0))
    do k=mink, maxk
      coff = [coff, ((-1)**k)*exp(cnst-fct(k+1)-fct(k+C+1)-fct(A-k+1)-fct(B-k+1))]
    end do
    
    ! makeup small d function
    db(:) = 0.0_real128
    do step=1,bsz
      cb =  cos(real(bet(step),kind=real64)/2)
      sb = -sin(real(bet(step),kind=real64)/2)
      kapp = 1
      do k=mink,maxk
        db(step) = db(step) + coff(kapp)*cb**(j-C-2*k)*sb**(C+2*k)
        kapp = kapp + 1
      end do
    end do
    
    ! makeup D function
    do concurrent (ia=1:asz, ib=1:bsz, ig=1:gsz)
      dwigner_d_(ia,ib,ig) = ea(ia)*eg(ig)*db(ib)
    end do

  end function dwigner_d_


  !>
  !> @returns 3x3 rotation matrix for 3D vector in Cartesian coordinate in real(real64) precision
  !> @note The matrix returned is the TRANSPOSE of the active rotation
  !>       R = Rz(alp) Ry(bet) Rz(gam), i.e. the passive rotation (rotation of the axes).
  !>       To rotate a vector v actively, use matmul(transpose(rot_cart(alp,bet,gam)), v).
  !>       rot_spin returns the active D(R) of the same R, so transpose(rot_cart) and rot_spin
  !>       rotate positions and spinors consistently: U^dag sigma_i U = sum_j R_ij sigma_j, U = rot_spin.
  !> @param alp alpha of Euler angles
  !> @param bet beta of Euler angles
  !> @param gam gamma of Euler angles
  pure function rot_cart(alp,bet,gam)
    implicit none
    ! returns
    real(real64) :: rot_cart(VDIM,VDIM)
    ! arguments 
    real(real64), intent(in) :: alp,bet,gam
    ! internal variables
    real(real64) :: c1,c2,c3,s1,s2,s3
    
    ! body
    c1 = cos(alp)
    c2 = cos(bet)
    c3 = cos(gam)
    s1 = sin(alp)
    s2 = sin(bet)
    s3 = sin(gam)
    
    rot_cart(1,1) = c1*c2*c3 - s1*s3
    rot_cart(1,2) = s1*c2*c3 + c1*s3
    rot_cart(1,3) = -s2*c3
    rot_cart(2,1) = -c1*c2*s3 - s1*c3
    rot_cart(2,2) = -s1*c2*s3 + c1*c3
    rot_cart(2,3) = s2*s3
    rot_cart(3,1) = c1*s2
    rot_cart(3,2) = s1*s2
    rot_cart(3,3) = c2
  end function rot_cart


  !> @returns 2x2 rotation matrix for 1/2 spinor in complex(real64) precision\
  !> @note the order of rows and columns are inverted from wigner_d (first row and column corresponds to spin up)
  !> @note The matrix returned is the active D^{1/2}(R) = exp(-i alp Sz) exp(-i bet Sy) exp(-i gam Sz),
  !>       not transposed.  The Cartesian counterpart is transpose(rot_cart(alp,bet,gam)).
  !> @param alp alpha of Euler angles
  !> @param bet beta of Euler angles
  !> @param gam gamma of Euler angles
  pure function rot_spin(alp,bet,gam)
    implicit none
    ! returns
    complex(real64) :: rot_spin(SDIM,SDIM)
    ! arguments
    real(real64), intent(in) :: alp
    real(real64), intent(in) :: bet
    real(real64), intent(in) :: gam
    ! local variables
    real(real64) :: c2,s2
    complex(real64) :: ea,eg

    ea = exp(cmplx(0.0d0,0.50d0,kind=real64)*alp)
    eg = exp(cmplx(0.0d0,0.50d0,kind=real64)*gam)
    c2 = cos(0.50d0*bet)
    s2 = sin(0.50d0*bet)
    rot_spin(SU,SU) = c2/(ea*eg)
    rot_spin(SU,SD) = -s2*eg/ea
    rot_spin(SD,SU) = s2*ea/eg
    rot_spin(SD,SD) = c2*ea*eg

  end function rot_spin



  !>
  !> @brief calculate Euler angles from rotation matrix for 3D vector in Cartesian representation
  !> @returns Euler angles alp, bet, gam corresponding to the given rotation matrix
  !> @param rot rotation matrix in Cartesian representation, in the convention of rot_cart
  !>            (the transpose of the active rotation); it inverts rot_cart
  subroutine euler_from_rot_cart(alp,bet,gam,rot)
    use, intrinsic :: ieee_arithmetic
    implicit none
    ! returns
    real(real64), intent(out) :: alp,bet,gam
    ! arguments
    real(real64), intent(in) :: rot(VDIM,VDIM)
    ! internal variables
    real(real64) :: dm

    ! body

    ! determine beta
    dm = rot(3,3)
    if(abs(dm) > 1.0_real64) then
      write(error_unit,'(A)') 'euler_from_rotation_cart: invalid value for Rx(3,3)'
      dm = sign(1.0_real64,dm)
    end if
    bet = acos(dm)

    ! determine alpha
    dm = rot(3,2)/rot(3,1)
    if(ieee_is_finite(dm)) then
      alp = atan(dm)
      if(rot(3,1) < 0.0_real64) alp = alp + PI
    else 
      write(error_unit,'(A)') 'euler_from_rotation_cart: invalid value for Rx(3,2)/Rx(3,1)'
      alp = -PI/2
      if(rot(3,2) > 0.0_real64) alp = alp + PI
    end if

    ! determine gamma
    dm = -rot(2,3)/rot(1,3)
    if(ieee_is_finite(dm)) then
      gam = atan(dm)
      if(rot(1,3) > 0.0_real64) gam = gam + PI
    else 
      write(error_unit,'(A)') 'euler_from_rotation_car: invalid value for Rx(3,2)/Rx(3,1)'
      gam = -PI/2
      if(rot(2,3) > 0.0_real64) gam = gam + PI
    end if

  end subroutine euler_from_rot_cart


  !>
  !> @brief calculate Clebsch-Gordan coefficients for given l1, l2, l3 and m3\
  !> @note l1, l2 and l3 are not doubled as in wigner's D function 
  !> @returns Clebsch-Gordan coefficients cg(m1,m2) := <l1,m1,l2,m2|l3,m3>
  !> @param l3 angular momentum of the total system
  !> @param m3 magnetic substitute of the total system
  !> @param l1 angular momentum of the first subsystem
  !> @param l2 angular momentum of the second subsystem
  pure function cg_coeff(l3,m3,l1,l2)
    implicit none
    ! returns
    real(real64), allocatable :: cg_coeff(:,:)
    ! arguments
    integer, intent(in) :: l3
    integer, intent(in) :: m3
    integer, intent(in) :: l1
    integer, intent(in) :: l2
    ! local variables
    integer :: jjcnst,max_a2,min_a2,a1,a2,jj
    real(real64) :: j1p,j1m,j2p,j2m,cnst
    
    ! body

    allocate(cg_coeff(l1+1,l2+1)) ! dimension is number of magnetic substitutes
    cg_coeff(:,:) = 0.0_real64
    jjcnst = l3*(l3+2) - l1*(l1+2) - l2*(l2+2)
    max_a2 = min(l2+1,(l1+l2+m3)/2 + 1)
    min_a2 = max(1,(-l1+l2+m3)/2 + 1)
    
    a2 = min_a2
    a1 = 2 + (l1+l2+m3)/2 - a2
    cg_coeff(a1,a2) = 1.0_real64
    if(a2+1 <= max_a2) then
      jj = (jjcnst - 2*(l1 - 2*(a1-1))*(l2 - 2*(a2-1)))/4
      j1m = sqrt(dble((a1-1)*(l1-a1+2)))
      j2p = sqrt(dble(a2*(l2-a2+1)))
      cg_coeff(a1-1,a2+1) = jj*cg_coeff(a1,a2)/(j1m*j2p)
    else
      return
    end if
    
    do a2=min_a2+1,max_a2-1
      a1 = 2 + (l1+l2+m3)/2 - a2
      jj = (jjcnst - 2*(l1 - 2*(a1-1))*(l2 - 2*(a2-1)))/4
      j1p = sqrt(dble(a1*(l1-a1+1)))
      j1m = sqrt(dble((a1-1)*(l1-a1+2)))
      j2p = sqrt(dble(a2*(l2-a2+1)))
      j2m = sqrt(dble((a2-1)*(l2-a2+2)))
      cg_coeff(a1-1,a2+1) = &
        & (jj*cg_coeff(a1,a2) - j1p*j2m*cg_coeff(a1+1,a2-1))/(j1m*j2p)
    end do
    cnst = sqrt(sum(cg_coeff(:,:)*cg_coeff(:,:)))
    cg_coeff(:,:) = cg_coeff(:,:)/cnst
  end function cg_coeff



  !>
  !> @brief calculate weights and abscissa of Gauss-Legendre quadrature\
  !> @note the weights and abscissa are automatically reallocated by this subroutine
  !> @returns weights and abscissa of Gauss-Legendre quadrature
  !> @param rank rank of the quadrature
  subroutine gauleg(rank,weight,abiscissa)
    implicit none
    ! arguments
    integer, intent(in) :: rank ! rank
    ! returns
    real(real64), allocatable, intent(out) :: weight(:) ! weight
    real(real64), allocatable, intent(out) :: abiscissa(:) ! abiscissa
    ! parameter
    real(real64), parameter :: EPS = 3.0d-14
    ! local variables
    integer :: m,i,j
    real(real64) :: p1,p2,p3,pp,z,z1

    ! body
    if(allocated(weight)) deallocate(weight)
    if(allocated(abiscissa)) deallocate(abiscissa)
    allocate(weight(rank),abiscissa(rank))
    
    ! calculate weight and abissias
    m = (rank+1)/2
    do i=1,m
      z = cos(PI*(i-0.25_real64)/(rank+0.50_real64))
      do while(.true.)
        p1 = 1.0_real64
        p2 = 0.0_real64
        do j=1,rank
          p3 = p2
          p2 = p1
          p1 = ((2*j - 1)*z*p2-(j-1)*p3)/j
        end do
        pp = rank*(z*p1 - p2)/(z*z - 1.0_real64)
        z1 = z
        z = z1 - p1/pp
        if(abs(z-z1) < EPS) exit
      end do
      abiscissa(i) = -z
      abiscissa(rank+1-i) = z
      weight(i) = 2.0_real64/((1.0_real64-z*z)*pp*pp)
      weight(rank+1-i) = weight(i)
    end do
  end subroutine gauleg
  

end module ModMathFunc
