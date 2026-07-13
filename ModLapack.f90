!>
!> @file ModLapack.f90
!> @brief Module for Lapack wrapper
!> @author Masaaki Kimura
!> @details
!>  matrix decomposition
!> - getrf: LU decomposition of a square matrix (real, complex)
!> - potrf: Cholesky decomposition of a positive definite matrix (real, complex)
!> - sytrf: Bunch-Kaufman decomposition of a symmetric matrix (real, complex)
!> - hetrf: Bunch-Kaufman decomposition of a hermite matrix (complex)
!> linear equation solvers
!> - getrs: solve linear equation by LU decomposed matrix (real, complex)
!> - potrs: solve linear equation by Cholesky decomposed matrix (real, complex)
!> - sytrs: solve linear equation by Bunch-Kaufman decomposed matrix (real, complex)
!> - hetrs: solve linear equation by Bunch-Kaufman decomposed matrix (complex)
!> matrix inversions
!> - getri: calculate inverse matrix from LU decomposed matrix (real, complex)
!> - potri: calculate inverse matrix from Cholesky decomposed matrix (real, complex)
!> - sytri: calculate inverse matrix from Bunch-Kaufman decomposed matrix (real, complex)
!> - hetri: calculate inverse matrix from Bunch-Kaufman decomposed matrix (complex)
!> - trtri: calculate inverse matrix of tri-diagonal matrix (real, complex)
!> matrix determinant
!> - getrd: calculate determinant from LU decomposed matrix (real, complex)
!> - potrd: calculate determinant from Cholesky decomposed matrix (real, complex)
!> eigenvalue problem
!> - syev: eigenvalue problem for real valued matrix
!> - heev: eigenvalue problem for complex valued matrix


module ModLapack
  use iso_fortran_env, only: real64, error_unit
  use ModFortIO, only: text_color
  implicit none
  private
  public :: Lapack

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! public module variables
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  character(len=1), public :: mlpk_error_color = 'r' ! set error text color

  ! Derived type which wraps Lapack subroutines
  type Lapack

    contains

    !////////////////////////////////////////////////////////////////////////////
    ! matrix decomposition
    !////////////////////////////////////////////////////////////////////////////
    ! LU decomposition of general matrix
    procedure, nopass, private :: dgetrf_, zgetrf_
    generic, public :: getrf => dgetrf_, zgetrf_ 
    ! Cholesky decomposition of a positive definite symmetric(hermite) matrix
    procedure, nopass, private :: dpotrf_, zpotrf_
    generic, public :: potrf => dpotrf_, zpotrf_
    ! Bunch-Kaufman decomposition of a symmetric(hermite) matrix
    procedure, nopass, private :: dsytrf_, zsytrf_
    generic, public :: sytrf => dsytrf_, zsytrf_
    procedure, nopass, public :: hetrf => zhetrf_

    !////////////////////////////////////////////////////////////////////////////
    ! linear equation solvers
    !////////////////////////////////////////////////////////////////////////////
    ! solve linear equation by LU decomposed matrix
    procedure, nopass, private :: dgetrs_, zgetrs_
    generic, public :: getrs => dgetrs_, zgetrs_
    ! solve linear equation by Cholesky decomposed matrix
    procedure, nopass, private :: dpotrs_, zpotrs_
    generic, public :: potrs => dpotrs_, zpotrs_
    ! solve linear equation by Bunch-Kaufman decomposed matrix
    procedure, nopass, private :: dsytrs_, zsytrs_
    generic, public :: sytrs => dsytrs_, zsytrs_
    procedure, nopass, public :: hetrs => zhetrs_

    !////////////////////////////////////////////////////////////////////////////
    ! matrix inversions
    !////////////////////////////////////////////////////////////////////////////
    ! calculate inverse matrix from LU decomposed matrix
    procedure, nopass, private :: dgetri_, zgetri_
    generic, public :: getri => dgetri_, zgetri_
    ! calculate inverse matrix from Cholesky decomposed matrix
    procedure, nopass, private :: dpotri_, zpotri_
    generic, public :: potri => dpotri_, zpotri_
    ! calculate inverse matrix from Bunch-Kaufman decomposed matrix
    procedure, nopass, private :: dsytri_, zsytri_
    generic, public :: sytri => dsytri_, zsytri_
    procedure, nopass, public :: hetri => zhetri_
    ! calculate inverse matrix of tridiagonal matrix
    procedure, nopass, private :: dtrtri_, ztrtri_
    generic, public :: trtri => dtrtri_, ztrtri_

    !////////////////////////////////////////////////////////////////////////////
    ! matrix determinant
    !////////////////////////////////////////////////////////////////////////////
    ! calculate determinant from LU decomposed matrix
    procedure, nopass, private :: dgetrd_, zgetrd_
    generic, public :: getrd => dgetrd_, zgetrd_
    procedure, nopass, private :: dpotrd_, zpotrd_
    generic, public :: potrd => dpotrd_, zpotrd_


    !////////////////////////////////////////////////////////////////////////////
    ! interfaces for eigenvalue problem
    !////////////////////////////////////////////////////////////////////////////
    ! eigenvalue problem for real valued matrix
    procedure, nopass, public :: syev => dsyev_
    ! eigenvalue problem for complex valued matrix
    procedure, nopass, public :: heev => zheev_

  end type Lapack
  
  contains  
 
  !////////////////////////////////////////////////////////////////////////////
  ! matrix decompositions
  !////////////////////////////////////////////////////////////////////////////

  !> @brief LU decompose general real valued matrix
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param piv pivot array
  logical function dgetrf_(src,piv)
    implicit none
    ! returns
    real(real64), intent(inout) :: src(:,:)
    integer, intent(out) :: piv(:)
    ! external
    external :: dgetrf
    ! local variables
    integer :: sz,inf
    
    ! body
    dgetrf_ = .false. ! default return value

    sz = size(src,dim=1)
    call dgetrf(sz,sz,src,sz,piv,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dgetrf_:')//' error in dgetrf'
      return
    end if
    dgetrf_ = .true.
  end function dgetrf_


  !> @brief LU decompose general complex valued matrix
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param piv pivot array
  logical function zgetrf_(src,piv)
    implicit none
    ! returns
    complex(real64), intent(inout) :: src(:,:)
    integer, intent(out) :: piv(:)
    ! external
    external :: zgetrf
    ! local variables
    integer :: sz, inf
    
    ! body
    zgetrf_ = .false. ! default return value

    sz = size(src,dim=1)
    call zgetrf(sz,sz,src,sz,piv,inf)
    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zgetrf_:')//' error in zgetrf'
      return
    end if
    zgetrf_ = .true.
  end function zgetrf_


  !> @brief choresky decomposition of real valuded
  !>            positive-definite symmetric matrix
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function dpotrf_(src,lower)
    implicit none
    ! returns
    real(real64), intent(inout) :: src(:,:)
    ! arguments
    logical, intent(in), optional :: lower
    ! external
    external :: dpotrf
    ! local variables
    integer :: sz,inf
    character(len=1) :: uplo
    
    ! body
    dpotrf_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if
    
    sz = size(src,dim=1)
    call dpotrf(uplo,sz,src,sz,inf)
    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dpotrf_:')//' error in dpotrf'
      return
    end if
    dpotrf_ = .true.
  end function dpotrf_
    

  
  !> @brief choresky decomposition for complex valuded
  !>            positive-definite hermite matrix
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function zpotrf_(src,lower)
    implicit none
    ! returns
    complex(real64), intent(inout) :: src(:,:)
    ! arguments
    logical, intent(in), optional :: lower
    ! external
    external :: zpotrf
    ! local variables
    integer :: sz,inf
    character(len=1) :: uplo
    
    ! body
    zpotrf_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(src,dim=1)
    call zpotrf(uplo,sz,src,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zpotrf_:')//' error in zpotrf'
      return
    end if
    zpotrf_ = .true.
  end function zpotrf_


  !> @brief Bunch-Kaufman decomposition of real valued symmetric matrix
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param piv pivot array with the same size as src
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function dsytrf_(src,piv,lower)
    implicit none
    ! returns
    real(real64), intent(inout) :: src(:,:)
    integer, intent(out) :: piv(:)
    ! arguments
    logical, intent(in), optional :: lower
    ! external
    external :: dsytrf
    ! local variables
    integer :: sz,inf
    real(real64) :: wk(size(src,dim=1))
    character(len=1) :: uplo

    ! body
    dsytrf_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(src,dim=1)
    call dsytrf(uplo,sz,src,sz,piv,wk,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dsytrf_:')//' error in dsytrf'
      return
    end if
    dsytrf_ = .true.
  end function dsytrf_


  !> @brief Bunch-Kaufman decomposition of complex valued symmetric matrix
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param piv pivot array with the same size as src
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function zsytrf_(src,piv,lower)
    implicit none
    ! returns
    complex(real64), intent(inout) :: src(:,:)
    integer, intent(out) :: piv(:)
    ! arguments
    logical, intent(in), optional :: lower
    ! external
    external :: zsytrf
    ! local variables
    integer :: sz,inf
    complex(real64) :: wk(size(src,dim=1))
    character(len=1) :: uplo

    ! body
    zsytrf_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(src,dim=1)
    call zsytrf(uplo,sz,src,sz,piv,wk,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zsytrf_:')//' error in zsytrf'
      return
    end if
    zsytrf_ = .true.
  end function zsytrf_


  !> @brief Bunch-Kaufman decomposition of complex valued hermite matrix
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param piv pivot array with the same size as src
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function zhetrf_(src,piv,lower)
    implicit none
    ! returns
    complex(real64), intent(inout) :: src(:,:)
    integer, intent(out) :: piv(:)
    ! arguments
    logical, intent(in), optional :: lower
    ! external
    external :: zhetrf
    ! local variables
    integer :: sz,inf
    complex(real64) :: wk(size(src,dim=1))
    character(len=1) :: uplo

    ! body
    zhetrf_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if
   
    sz = size(src,dim=1)
    call zhetrf(uplo,sz,src,sz,piv,wk,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zhetrf_:')//' error in zhetrf'
      return
    end if
    zhetrf_ = .true.
  end function zhetrf_


  !////////////////////////////////////////////////////////////////////////////
  ! linear equation solvers
  !////////////////////////////////////////////////////////////////////////////

  !> @brief solve linear equation from lu of real valued general matrix
  !> @return true if successful, false otherwise
  !> @param a lu decomposed matrix
  !> @param piv pivot array from lu decomposition
  !> @param b right-hand side matrix
  !> @param trans if 'T' or 'C', solve transposed problem
  logical function dgetrs_(a,piv,b,trans)
    implicit none
    ! returns
    real(real64), intent(inout) :: b(:,:)
    ! arguments
    real(real64), intent(in) :: a(:,:)
    integer, intent(in) :: piv(:)
    character(len=*), intent(in), optional :: trans
    ! external
    external :: dgetrs
    ! local variables
    integer :: sz1,sz2,inf
    
    ! body
    dgetrs_ = .false. ! default return value

    sz1 = size(a,dim=1)
    sz2 = size(b,dim=2)
    if(present(trans)) then
      if(trans/='N' .and. trans/='T' .and. trans/= 'C') then
        print*, 'lapack%dgetrs_: invalid value for trans'
        return
      end if
      call dgetrs(trans,sz1,sz2,a,sz1,piv,b,sz1,inf)
    else
      call dgetrs('N',sz1,sz2,a,sz1,piv,b,sz1,inf)
    end if

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dgetrs_:')//' error in dgetrs'
      return
    end if
    dgetrs_ = .true.
  end function dgetrs_



  !> @brief Solves linear equation from lu of complex valued general matrix
  !> @return true if successful, false otherwise
  !> @param a lu decomposed matrix
  !> @param piv pivot array from lu decomposition
  !> @param b right-hand side matrix
  !> @param trans if 'T' or 'C', solve transposed problem
  logical function zgetrs_(a,piv,b,trans)
    implicit none
    ! returns
    complex(real64), intent(inout) :: b(:,:)
    ! arguments
    complex(real64), intent(in) :: a(:,:)
    integer, intent(in) :: piv(:)
    character(len=*), intent(in), optional :: trans
    ! external
    external :: zgetrs    
    ! local variables
    integer :: sz1,sz2,inf
    
    ! body
    zgetrs_ = .false. ! default return value

    sz1 = size(a,dim=1)
    sz2 = size(b,dim=2)
    if(present(trans)) then
      if(trans/= 'N' .and. trans/= 'T' .and. trans/= 'C') then
        print*, 'lapack%zgetrs_: invalid value for trans'
        return
      end if
      call zgetrs(trans,sz1,sz2,a,sz1,piv,b,sz1,inf)
    else
      call zgetrs('N',sz1,sz2,a,sz1,piv,b,sz1,inf)
    end if

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zgetrs_:')//' error in zgetrs'
      return
    end if
    zgetrs_ = .true.
  end function zgetrs_
  
  
  !> @brief Solve linear equation from choresky decomposition of positive-definite
  !>        real-valued symmetric matrix
  !> @return true if successful, false otherwise
  !> @param a choresky decomposed matrix
  !> @param b right-hand side matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function dpotrs_(a,b,lower)
    implicit none
    ! returns
    real(real64), intent(inout) :: b(:,:)
    ! argument
    real(real64), intent(in) :: a(:,:)
    logical, intent(in), optional :: lower
    ! external
    external :: dpotrs
    ! local variables
    integer :: sz,inf
    character(len=1) :: uplo
    
    ! body
    dpotrs_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(a,dim=1)
    call dpotrs(uplo,sz,sz,a,sz,b,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dpotrs_:')//' error in dpotrs'
      return
    end if
    dpotrs_ = .true.
  end function dpotrs_


  !> @brief Solves linear equation from choresky decomposition of positive-definite
  !>        complex-valued hermite matrix
  !> @return true if successful, false otherwise
  !> @param a choresky decomposed matrix
  !> @param b right-hand side matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function zpotrs_(a,b,lower)
    implicit none
    ! returns
    complex(real64), intent(inout) :: b(:,:)
    ! argument
    complex(real64), intent(in) :: a(:,:)
    logical, intent(in), optional :: lower
    ! external
    external :: zpotrs
    ! local variables
    integer :: sz,inf
    character(len=1) :: uplo
    
    ! body
    zpotrs_ = .false. ! default return value
    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(a,dim=1)
    call zpotrs(uplo,sz,sz,a,sz,b,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zpotrs_:')//' error in zpotrs'
      return
    end if
    zpotrs_ = .true.
  end function zpotrs_
  
  

  !> @brief Solves lienar equation from Bunch-Kaufman decomposed real valued symmetric matrix
  !> @return true if successful, false otherwise
  !> @param a Bunch-Kaufman decomposed matrix
  !> @param piv pivot array from Bunch-Kaufman decomposition
  !> @param b right-hand side matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function dsytrs_(a,piv,b,lower)
    implicit none
    ! return
    real(real64), intent(inout) :: b(:,:)
    ! argument
    real(real64), intent(in) :: a(:,:)
    integer, intent(in) :: piv(:)
    logical, intent(in), optional :: lower
    ! external
    external :: dsytrs    
    ! local variables
    integer :: sz,inf
    character(len=1) :: uplo
    
    ! body
    dsytrs_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(a,dim=1)
    call dsytrs(uplo,sz,sz,a,sz,piv,b,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dsytrs_:')//' error in dsytrs'
      return
    end if
    dsytrs_ = .true.
  end function dsytrs_


  !> @brief Solves lienar equation from Bunch-Kaufman decomposed complex valued symmetric matrix
  !> @return true if successful, false otherwise
  !> @param a Bunch-Kaufman decomposed matrix
  !> @param piv pivot array from Bunch-Kaufman decomposition
  !> @param b right-hand side matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function zsytrs_(a,piv,b,lower)
    implicit none
    ! return
    complex(real64), intent(inout) :: b(:,:)
    ! argument
    complex(real64), intent(in) :: a(:,:)
    integer, intent(in) :: piv(:)
    logical, intent(in), optional :: lower
    ! external
    external :: zsytrs
    ! local variables
    integer :: sz,inf
    character(len=1) :: uplo
    
    ! body
    zsytrs_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(a,dim=1)
    call zsytrs(uplo,sz,sz,a,sz,piv,b,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zsytrs_:')//' error in zsytrs'
      return
    end if
    zsytrs_ = .true.
  end function zsytrs_


  !> @brief Solves lienar equation from Bunch-Kaufman decomposed complex valued hermite matrix
  !> @return true if successful, false otherwise
  !> @param a Bunch-Kaufman decomposed matrix
  !> @param piv pivot array from Bunch-Kaufman decomposition
  !> @param b right-hand side matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function zhetrs_(a,piv,b,lower)
    implicit none
    ! return
    complex(real64), intent(inout) :: b(:,:)
    ! argument
    complex(real64), intent(in) :: a(:,:)
    integer, intent(in) :: piv(:)
    logical, intent(in), optional :: lower
    ! external
    external :: zhetrs    
    ! local variables
    integer :: sz,inf
    character(len=1) :: uplo
    
    ! body
    zhetrs_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(a,dim=1)
    call zhetrs(uplo,sz,sz,a,sz,piv,b,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zhetrs_:')//' error in zhetrs'
      return
    end if
    zhetrs_ = .true.
  end function zhetrs_

  
  
  !////////////////////////////////////////////////////////////////////////////
  ! matrix inversions
  !////////////////////////////////////////////////////////////////////////////
  
  !> @brief Calculates inverse matrix of real valued general matrix from LU decomposed matrix
  !> @return true if successful, false otherwise
  !> @param src LU decomposed matrix
  !> @param piv pivot array from LU decomposition
  logical function dgetri_(src,piv)
    implicit none
    ! returns
    real(real64), intent(inout) :: src(:,:)
    ! arguments
    integer, intent(in) :: piv(:)
    ! external
    external :: dgetri
    ! local variables
    integer :: sz,inf
    real(real64) :: wk(size(src,dim=1))
    
    ! body
    dgetri_ = .false. ! default return value

    sz = size(src,dim=1)
    call dgetri(sz,src,sz,piv,wk,sz,inf)    

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dgetri_:')//' error in dgetri'
      return
    end if
    dgetri_ = .true.
  end function dgetri_


  !> @brief Calculates inverse matrix of complex valued general matrix from LU decomposed matrix
  !> @return true if successful, false otherwise
  !> @param src LU decomposed matrix
  !> @param piv pivot array from LU decomposition
  logical function zgetri_(src,piv)
    implicit none
    ! returns
    complex(real64), intent(inout) :: src(:,:)
    ! arguments
    integer, intent(in) :: piv(:)
    ! external
    external :: zgetri    
    ! local variables
    integer :: sz,inf
    complex(real64) :: wk(size(src,dim=1))
    
    ! body
    zgetri_ = .false. ! default return value

    sz = size(src,dim=1)
    call zgetri(sz,src,sz,piv,wk,sz,inf)    

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zgetri_:')//' error in zgetri'
      return
    end if
    zgetri_ = .true.
  end function zgetri_


  !> @brief Calculates inverse matrix of real valued positive definite symmetric matrix 
  !>        from choresky decomposed matrix
  !> @return true if successful, false otherwise
  !> @param src Cholesky decomposed matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function dpotri_(src,lower)
    implicit none
    ! returns
    real(real64), intent(inout) :: src(:,:)
    ! arguments
    logical, intent(in), optional :: lower    
    ! external
    external :: dpotri    
    ! local variables
    integer :: sz,i,j,inf
    character(len=1) :: uplo
    
    ! body
    dpotri_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(src,dim=1)
    call dpotri(uplo,sz,src,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dpotri_:')//' error in dpotri'
      return
    end if

    if(uplo=='U') then
      ! copy upper triangle to lower triangle
      do i=1,sz-1
        do j=i+1,sz
          src(j,i) = src(i,j)
        end do
      end do
    else
      ! copy lower triangle to upper triangle
      do i=1,sz-1
        do j=i+1,sz
          src(i,j) = src(j,i)
        end do
      end do
    end if
    dpotri_ = .true.
  end function dpotri_


  !> @brief Calculates inverse matrix of complex valued positive definite hermite matrix
  !>        from choresky decomposed matrix
  !> @return true if successful, false otherwise
  !> @param src Cholesky decomposed matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle

  logical function zpotri_(src,lower)
    implicit none
    ! returns
    complex(real64), intent(inout) :: src(:,:)
    ! arguments
    logical, intent(in), optional :: lower
    ! external
    external :: zpotri    
    ! local variables
    integer :: sz,i,j,inf
    character(len=1) :: uplo
    
    ! body
    zpotri_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(src,dim=1)
    call zpotri(uplo,sz,src,sz,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zpotri_:')//' error in zpotri'
      return
    end if

    if(uplo=='U') then
      ! copy upper triangle to lower triangle
      do i=1,sz-1
        do j=i+1,sz
          src(j,i) = conjg(src(i,j))
        end do
      end do
    else
      ! copy lower triangle to upper triangle
      do i=1,sz-1
        do j=i+1,sz
          src(i,j) = conjg(src(j,i))
        end do
      end do
    end if
    zpotri_ = .true.
  end function zpotri_

  
  !> @brief Calculates inverse matrix of real vlaued symmetric matrix by Bunch-Kaufman decomposition
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param piv pivot array from Bunch-Kaufman decomposition
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function dsytri_(src,piv,lower)
    implicit none
    ! returns
    real(real64), intent(inout) :: src(:,:)
    ! arguments
    integer, intent(in) :: piv(:)
    logical, intent(in), optional :: lower
    ! external
    external :: dsytri    
    ! local variables
    integer :: sz,inf
    real(real64) :: wk(2*size(src,dim=1))
    character(len=1) :: uplo
    
    ! body
    dsytri_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(src,dim=1)
    call dsytri(uplo,sz,src,sz,piv,wk,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dsytri_:')//' error in dsytri'
      return
    end if
    dsytri_ = .true.
  end function dsytri_


  !> @brief Calculate inverse matrix of complex vlaued symmetric matrix 
  !>        by Bunch-Kaufman decomposition
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param piv pivot array from Bunch-Kaufman decomposition
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function zsytri_(src,piv,lower)
    implicit none
    ! returns
    complex(real64), intent(inout) :: src(:,:)
    ! arguments
    integer, intent(in) :: piv(:)
    logical, intent(in), optional :: lower
    ! external
    external :: zsytri    
    ! local variables
    integer :: sz,inf
    complex(real64) :: wk(2*size(src,dim=1))
    character(len=1) :: uplo
    
    ! body
    zsytri_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(src,dim=1)
    call zsytri(uplo,sz,src,sz,piv,wk,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zsytri_:')//' error in zsytri'
      return
    end if
    zsytri_ = .true.
  end function zsytri_



  !> @brief Calculates inverse matrix of complex vlaued hermite matrix
  !>        by Bunch-Kaufman decomposition
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param piv pivot array from Bunch-Kaufman decomposition
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function zhetri_(src,piv,lower)
    implicit none
    ! returns
    complex(real64), intent(inout) :: src(:,:)
    ! arguments
    integer, intent(in) :: piv(:)
    logical, intent(in), optional :: lower
    ! external
    external :: zhetri    
    ! local variables
    integer :: sz,inf
    complex(real64) :: wk(2*size(src,dim=1))
    character(len=1) :: uplo
    
    ! body
    zhetri_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    sz = size(src,dim=1)
    call zhetri(uplo,sz,src,sz,piv,wk,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zhetri_:')//' error in zhetri'
      return
    end if
    zhetri_ = .true.
  end function zhetri_


  !> @brief Calculates inverse matrix of real valued triangular matrix
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  !> @param unity if true, unit diagonal is assumed
  logical function dtrtri_(src,lower,unity)
    implicit none
    ! returns
    real(real64), intent(inout) :: src(:,:)
    ! arguments
    logical, intent(in), optional :: lower
    logical, intent(in), optional :: unity
    ! external
    external :: dtrtri
    ! local variables
    integer :: sz,inf
    character(len=1) :: uplo,diag
    

    ! body
    dtrtri_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    diag = 'N'  ! non-unit diagonal (default)
    if(present(unity)) then
      if(unity) diag = 'U'  ! unit diagonal
    end if

    sz = size(src,dim=1)
    call dtrtri(uplo, diag, sz, src, sz, inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dtrtri_:')//' error in dtrtri'
      return
    end if
    dtrtri_ = .true.
  end function dtrtri_

  !> @brief Calculates inverse matrix of complex valued triangular matrix
  !> @return true if successful, false otherwise
  !> @param src input matrix
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  !> @param unity if true, unit diagonal is assumed
  logical function ztrtri_(src,lower,unity)
    implicit none
    ! returns
    complex(real64), intent(inout) :: src(:,:)
    ! arguments
    logical, intent(in), optional :: lower
    logical, intent(in), optional :: unity
    ! external
    external :: ztrtri
    ! local variables
    integer :: sz,inf
    character(len=1) :: uplo,diag
    

    ! body
    ztrtri_ = .false. ! default return value

    uplo = 'U'  ! upper triangle is decomposed (default)
    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    end if

    diag = 'N'  ! non-unit diagonal (default)
    if(present(unity)) then
      if(unity) diag = 'U'  ! unit diagonal
    end if

    sz = size(src,dim=1)
    call ztrtri(uplo, diag, sz, src, sz, inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%ztrtri_:')//' error in ztrtri'
      return
    end if
    ztrtri_ = .true.
  end function ztrtri_



  
  !////////////////////////////////////////////////////////////////////////////
  ! determinants
  !////////////////////////////////////////////////////////////////////////////
  
  !> @brief Calculates determinant of real valued general matrix by LU decomposition
  !> @return determinant
  !> @param src input matrix
  !> @param piv pivot array from LU decomposition
  pure real(real64) function dgetrd_(src,piv)
    implicit none
    ! arguments
    real(real64), intent(in) :: src(:,:)
    integer, intent(in) :: piv(:)
    ! local variables
    integer :: sz,i

    ! body
    sz = size(src,dim=1)
    dgetrd_ = 1.0_real64
    do i=1,sz
      dgetrd_ = dgetrd_*src(i,i)
      if(piv(i) /= i) dgetrd_ = -1*dgetrd_
    end do
  end function dgetrd_


  !> @brief Calculates determinant of complex valued general matrix by LU decomposition
  !> @return determinant
  !> @param src input matrix
  !> @param piv pivot array from LU decomposition
  pure complex(real64) function zgetrd_(src,piv)
    implicit none
    ! arguments
    complex(real64), intent(in) :: src(:,:)
    integer, intent(in) :: piv(:)
    ! local variables
    integer :: sz,i

    ! body
    sz = size(src,dim=1)
    zgetrd_ = 1.0_real64
    do i=1,sz
      zgetrd_ = zgetrd_*src(i,i)
      if(piv(i) /= i) zgetrd_ = -1*zgetrd_
    end do
  end function zgetrd_


  !> @brief Calculates determinant of real valued positive definite
  !>               symmetric matrix by Cholesky decomposition
  !> @return determinant
  !> @param src Cholesky decomposed matrix
  pure real(real64) function dpotrd_(src)
    implicit none
    ! arguments
    real(real64), intent(in) :: src(:,:)
    ! local variables
    integer :: sz,i
    
    ! body
    sz = size(src,dim=1)
    dpotrd_ = 1.0_real64
    do i=1,sz
      dpotrd_ = dpotrd_*src(i,i)
    end do
    dpotrd_ = dpotrd_*dpotrd_
  end function dpotrd_


  !> @brief Calculates determinant of complex valued positive definite
  !>               hermite matrix by Cholesky decomposition
  !> @return determinant
  !> @param src Cholesky decomposed matrix
  pure complex(real64) function zpotrd_(src)
    implicit none
    ! arguments
    complex(real64), intent(in) :: src(:,:)
    ! local variables
    integer :: sz,i
    
    ! body
    sz = size(src,dim=1)
    zpotrd_ = 1.0_real64
    do i=1,sz
      zpotrd_ = zpotrd_*src(i,i)
    end do
    zpotrd_ = zpotrd_*zpotrd_
  end function zpotrd_
  
  

  !////////////////////////////////////////////////////////////////////////////
  ! eigenvalue problem
  !////////////////////////////////////////////////////////////////////////////

  !> @brief Solves eivenvalue problem for real valued symmetric matrix
  !> @return true if successful, false otherwise
  !> @param a input matrix
  !> @param w eigenvalues
  !> @param jobz if 'N', only eigenvalues are calculated
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function dsyev_(a,w,jobz,lower)
    implicit none
    ! return
    real(real64), intent(inout) :: a(:,:)
    real(real64), intent(out) :: w(:)
    ! argument
    character, intent(in) :: jobz
    logical, intent(in), optional :: lower
    ! externals
    external :: dsyev
    ! local variables
    integer :: sz,lwk,inf
    real(real64) :: wk(max(3*size(a,dim=1)-1,1))
    character(len=1) :: uplo
    
    ! body
    dsyev_ = .false. ! default return value

    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    else
      uplo = 'U'  ! upper triangle is decomposed (default)
    end if

    sz = size(a,dim=1)
    lwk = max(3*sz -1,1)
    call dsyev(jobz,uplo,sz,a,sz,w,wk,lwk,inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%dsyev_:')//' error in dsyev'
      return
    end if
    dsyev_ = .true.
  end function dsyev_


  !> @brief Solves eigenvalue problem for complex valued hermite matrix
  !> @return true if successful, false otherwise
  !> @param a input matrix
  !> @param w eigenvalues
  !> @param jobz if 'N', only eigenvalues are calculated
  !> @param lower if true, decomposes lower triangle, otherwise upper triangle
  logical function zheev_(a,w,jobz,lower)
    implicit none
    ! return
    complex(real64), intent(inout) :: a(:,:)
    real(real64), intent(out) :: w(:)
    ! argument
    character, intent(in) :: jobz
    logical, intent(in), optional :: lower
    ! externals
    external :: zheev
    ! local variables
    integer sz,lwk,inf
    real(real64) :: rwk(max(1,3*size(a,dim=1)-2))
    complex(real64) :: wk(max(1,2*size(a,dim=1)-1))
    character(len=1) :: uplo
    
    ! body
    zheev_ = .false. ! default return value

    if(present(lower)) then
      if(lower) uplo = 'L'  ! lower triangle is decomposed
    else
      uplo = 'U'  ! upper triangle is decomposed (default)
    end if

    sz = size(a,dim=1)
    lwk = 2*sz-1
    call zheev(jobz,uplo,sz,a,sz,w,wk,lwk,rwk, inf)

    if(inf/=0) then
      write(error_unit,'(A)') text_color(mlpk_error_color,'lapack%zheev_:')//' error in zheev'
      return
    end if
    zheev_ = .true.
  end function zheev_

end module ModLapack

