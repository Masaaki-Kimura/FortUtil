!> @file ModSpinParity.f90
!! @brief Module for spin-parity symbol
!! @author Masaaki Kimura
module ModSpinParity
  use iso_fortran_env, only: real64, error_unit
  use ModNuclConst, only: POS, NEG, PARITY_CHAR
  use ModFortIO, only: MAXBUF, text_color
  implicit none
  private
  ! class for spin parity
  public :: SpinParity 
  ! logical operators
  public :: operator(==)
  public :: operator(/=)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! public module variables
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  character(len=1), public :: msp_error_color = 'r' ! set error text color

  !> @brief stores the information of spin and parity
  type SpinParity
    ! initialize with abnormal values
    integer :: j = -1     ! angular momentum
    integer :: dim = -1   ! dimension of magnetic substate
    integer :: p =-1      ! parity (1:POS, 2:NEG)
    
  contains
    procedure, public :: m => sp_m  ! return magnetic substate for given index
    procedure, public :: m2index => sp_m2index ! return index for given magnetic substate
    procedure, public :: str => sp_str ! return strings of the spin-parity
    
    ! check if the spin-parity satisfies the triangle condition
    procedure, public :: is_triangle => sp_is_triangle
    ! check if the spin-parity is natural parity state
    procedure, public :: is_natural => sp_is_natural
    
    ! return index for positive parity
    procedure, nopass, public :: pos => sp_pos
    ! return index for negative parity
    procedure, nopass, public :: neg => sp_neg
    ! return zero angular momentum
    procedure, nopass, public :: zero => sp_zero

  end type SpinParity
  
  ! interface for SpinParity constructor
  interface SpinParity
    module procedure sp_init_by_int, sp_init_by_str
  end interface SpinParity


  
  ! interface for SpinParity operator(==)
  interface operator(==)
    module procedure sp_is_equal
  end interface operator(==)

  ! interface for SpinParity operator(/=)
  interface operator(/=)
    module procedure sp_is_unequal
  end interface operator(/=)

  
  contains

  !> @brief Sets spin-parity from string input
  !! @return SpinParity object
  !! @param src Input string
  type(SpinParity) function sp_init_by_str(src)
    implicit none
    ! argument
    character(len=*) :: src
    ! local variables
    integer :: is
    character(:), allocatable :: buf

    ! body
    ! check the length of input
    buf = trim(adjustl(src))
    if (len(buf) < 2) then
      write(error_unit,*) text_color(msp_error_color,'sp_init_by_str:')//' Error in input string, ',buf
      return
    end if
    
    ! read angular momentum
    read(buf(1:len(buf)-1),*,iostat=is) sp_init_by_str%j
    if(is /= 0 .or. sp_init_by_str%j < 0) then
      write(error_unit,*) text_color(msp_error_color,'sp_init_by_str:')//' Error in input string, ',buf
      return
    end if
    
    ! set dimension of magnetic substate
    sp_init_by_str%dim = sp_init_by_str%j + 1
    
    ! read parity
    select case(buf(len(buf):len(buf)))
      case ('+')
        sp_init_by_str%p = POS
      case ('-')
        sp_init_by_str%p = NEG
      case default
        write(error_unit,*) text_color(msp_error_color,'sp_init_by_str:')//' Error in input string, ',trim(buf)
        return
    end select
  end function sp_init_by_str

  
  !> @brief Sets spin-parity from integer input
  !! @return SpinParity object
  !! @param j Angular momentum
  !! @param p Parity (1:POS, 2:NEG)
  type(SpinParity) function sp_init_by_int(j,p)
    implicit none
    ! argument
    integer, intent(in) :: j
    integer, intent(in) :: p

    ! body
    if(j<0 .or. (p/=POS .and. p/=NEG)) then
      write(error_unit,'(A)') text_color(msp_error_color,'sp_init_by_str:')//' invalid input for j or p'
      return
    end if
    ! set angular momentum and parity
    sp_init_by_int%j = j;   sp_init_by_int%p = p
    
    ! set dimension of magnetic substate
    sp_init_by_int%dim = sp_init_by_int%j + 1
  end function sp_init_by_int

  
  !> @brief Checks if two SpinParity objects are equal
  !! @return .true. if lhs and rhs are equal, .false. otherwise
  !! @param lhs Left-hand side SpinParity object
  !! @param rhs Right-hand side SpinParity object
  elemental logical function sp_is_equal(lhs,rhs)
    implicit none
    ! arguments
    type(SpinParity), intent(in) :: lhs
    type(SpinParity), intent(in) :: rhs
    ! body
    sp_is_equal = lhs%j == rhs%j .and. lhs%p == rhs%p
  end function sp_is_equal


  !> @brief Checks if lhs and rhs spin-parity are unequal
  !! @return .true. if lhs and rhs are unequal, .false. otherwise
  !! @param lhs Left-hand side SpinParity object
  !! @param rhs Right-hand side SpinParity object
  elemental logical function sp_is_unequal(lhs,rhs)
    implicit none
    ! arguments
    type(SpinParity), intent(in) :: lhs
    type(SpinParity), intent(in) :: rhs
    ! body
    sp_is_unequal = .not.(sp_is_equal(lhs,rhs))
  end function sp_is_unequal

  
  !> @brief  Returns string representation of the spin-parity
  !! @return String representation of the spin-parity
  pure function sp_str(this)
    implicit none
    ! returns
    character(:), allocatable :: sp_str
    ! arguments
    class(SpinParity), intent(in) :: this
    ! local variables
    character(len=MAXBUF) :: jstr

    ! write j and parity
    write(jstr ,'(I0)') this%j
    sp_str = trim(adjustl(jstr))//PARITY_CHAR(this%p)  
  end function sp_str
  
  
  !> @brief Returns magnetic substate from index\
  !! @note There is no check for the index range (1 <= i <= dim)
  !! @return Magnetic substate
  !! @param i Index
  elemental integer function sp_m(this,i)
    implicit none
    class(SpinParity), intent(in) :: this
    ! argument
    integer, intent(in) :: i
    ! body
    sp_m = -this%j + 2*(i-1)
  end function sp_m

  !> @brief Returns index for magnetic substate\
  !! @note There is no check for the magnetic substate range (-j <= m <= j)  
  !! @return Index corresponding to the given magnetic substate
  !! @param m Magnetic substate
  elemental integer function sp_m2index(this,m)
    implicit none
    class(SpinParity), intent(in) :: this
    ! argument
    integer, intent(in) :: m
    ! body
    sp_m2index = (m+this%j)/2+1
  end function sp_m2index


  !> @brief Checks if the spin-parity combination satisfies the triangle condition
  !! @return .true. if the triangle condition is satisfied, .false. otherwise
  !! @param jp2 Second SpinParity object
  !! @param jp3 Third SpinParity object
  elemental logical function sp_is_triangle(this,jp2,jp3)
    implicit none
    ! arguments
    class(SpinParity), intent(in) :: this ! spin-parity
    type(SpinParity), intent(in) :: jp2
    type(SpinParity), intent(in) :: jp3
    ! local variables
    integer :: j1, j2, j3
    ! body
    j1 = this%j;   j2 = jp2%j;    j3 = jp3%j
    sp_is_triangle = (mod(this%p + jp2%p + jp3%p, 2) == 1 ) .and. & ! parity
                     (abs(j1-j2) <= j3 .and. j3 <= j1 + j2) .and. & ! triangle
                     (mod(j1 + j2 + j3, 2) == 0 ) ! even sum of (2j)'s 
  end function sp_is_triangle


  !> @brief Checks if the spin-parity is natural parity state
  !! @return .true. if the spin-parity is natural parity state, .false. otherwise
  elemental logical function sp_is_natural(this)
    implicit none
    ! arguments
    class(SpinParity), intent(in) :: this
    ! local variables
    integer :: m
    ! body
    ! check if this is natural parity state
    m = mod(this%j,4)
    sp_is_natural = (m==0 .and. this%p==POS) .or. (m==2 .and. this%p==NEG )
  end function sp_is_natural


  !> @brief Returns index for positive parity
  !! @return Index for positive parity
  elemental integer function sp_pos()
    implicit none
    ! returns
    sp_pos = POS
  end function sp_pos

  !> @brief Returns index for negative parity
  !! @return Index for negative parity
  elemental integer function sp_neg()
    implicit none
    ! returns
    sp_neg = NEG
  end function sp_neg

  !> @brief Returns zero angular momentum
  !! @return Zero 
  elemental integer function sp_zero()
    implicit none
    sp_zero = 0
  end function sp_zero
    

  
end module ModSpinParity
