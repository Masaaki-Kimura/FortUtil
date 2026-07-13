program TestSpecialFunctions
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use mod_special_functions
  implicit none

  type(special_functions) :: sf

  call sf%timestamp()
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPECIAL_FUNCTIONS_TEST:'
  write ( *, '(a)' ) '  Adapted for mod_special_functions'
  write ( *, '(a)' ) '  Test selected SPECIAL_FUNCTIONS routines.'

  call airya_test(sf)
  call beta_test(sf)
  call cisia_test(sf)
  call cisib_test(sf)
  call cjy01_test(sf)
  call comelp_test(sf)
  call hygfx_test(sf)
  call sphj_test(sf)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPECIAL_FUNCTIONS_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call sf%timestamp()

contains

  subroutine airya_test(sf)
    type(special_functions), intent(in) :: sf

    integer(int32), parameter :: test_num = 4
    real(real64) :: ad, ai, bd, bi, x
    real(real64), dimension(test_num) :: x_test = (/ &
      0.0d0, 10.0d0, 20.0d0, 30.0d0 /)
    integer(int32) :: i

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AIRYA_TEST'
    write ( *, '(a)' ) '  Test AIRYA'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    x      Ai(x)         Bi(x)         Ai''(x)        Bi''(x)'
    write ( *, '(a)' ) ' '

    do i = 1, test_num
      x = x_test(i)
      call sf%airya(x, ai, bi, ad, bd)
      write ( *, '(1x,f5.1,4g16.8)' ) x, ai, bi, ad, bd
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    x     Ai(-x)        Bi(-x)        Ai''(-x)       Bi''(-x)'
    write ( *, '(a)' ) ' '

    do i = 1, test_num
      x = x_test(i)
      call sf%airya(-x, ai, bi, ad, bd)
      write ( *, '(1x,f5.1,4g16.8)' ) x, ai, bi, ad, bd
    end do
  end subroutine airya_test

  subroutine beta_test(sf)
    type(special_functions), intent(in) :: sf

    integer(int32), parameter :: test_num = 3
    real(real64) :: bt, p, q
    real(real64), dimension(test_num) :: p_test = (/ 1.5d0, 2.5d0, 1.5d0 /)
    real(real64), dimension(test_num) :: q_test = (/ 2.0d0, 2.0d0, 3.0d0 /)
    integer(int32) :: test

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_TEST:'
    write ( *, '(a)' ) '  Test BETA.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    p       q           B(p,q)'
    write ( *, '(a)' ) '  ---------------------------------'

    do test = 1, test_num
      p = p_test(test)
      q = q_test(test)
      call sf%beta(p, q, bt)
      write ( *, '(2x,f5.1,3x,f5.1,d20.10)' ) p, q, bt
    end do
  end subroutine beta_test

  subroutine cisia_test(sf)
    type(special_functions), intent(in) :: sf

    integer(int32), parameter :: test_num = 6
    real(real64) :: ci, si, x
    real(real64), dimension(test_num) :: x_test = (/ &
      0.0d0, 5.0d0, 10.0d0, 20.0d0, 30.0d0, 40.0d0 /)
    integer(int32) :: test

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CISIA_TEST'
    write ( *, '(a)' ) '  CISIA computes the cosine and sine integrals.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   x        ci(x)           si(x)'
    write ( *, '(a)' ) '------------------------------------'

    do test = 1, test_num
      x = x_test(test)
      call sf%cisia(x, ci, si)
      write ( *, '(1x,f5.1,g16.8,g16.8)' ) x, ci, si
    end do
  end subroutine cisia_test

  subroutine cisib_test(sf)
    type(special_functions), intent(in) :: sf

    integer(int32), parameter :: test_num = 6
    real(real64) :: ci, si, x
    real(real64), dimension(test_num) :: x_test = (/ &
      0.0d0, 5.0d0, 10.0d0, 20.0d0, 30.0d0, 40.0d0 /)
    integer(int32) :: test

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CISIB_TEST'
    write ( *, '(a)' ) '  CISIB computes the cosine and sine integrals.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   x        ci(x)           si(x)'
    write ( *, '(a)' ) '------------------------------------'

    do test = 1, test_num
      x = x_test(test)
      call sf%cisib(x, ci, si)
      write ( *, '(1x,f5.1,g16.8,g16.8)' ) x, ci, si
    end do
  end subroutine cisib_test

  subroutine cjy01_test(sf)
    type(special_functions), intent(in) :: sf

    complex(real64) :: cbj0, cbj1, cby0, cby1
    complex(real64) :: cdj0, cdj1, cdy0, cdy1, z
    real(real64) :: x, y

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CJY01_TEST'
    write ( *, '(a)' ) '  Test CJY01'
    write ( *, '(a)' ) ' '

    x = 4.0d0
    y = 2.0d0
    z = cmplx(x, y, kind=real64)
    write ( *, '(a,g16.8,a,g16.8)' ) '  Z = ', x, ' + i * ', y

    call sf%cjy01(z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1)

    write ( *, * )
    write ( *, * ) '  n      Re[Jn(z)]       Im[Jn(z)]       Re[Jn''(z)]      Im[Jn''(z)]'
    write ( *, * ) ' --------------------------------------------------------------------'
    write ( *, '(6x,4d16.8)' ) cbj0, cdj0
    write ( *, '(6x,4d16.8)' ) cbj1, cdj1

    write ( *, * )
    write ( *, * ) '  n      Re[Yn(z)]       Im[Yn(z)]       Re[Yn''(z)]      Im[Yn''(z)]'
    write ( *, * ) ' --------------------------------------------------------------------'
    write ( *, '(6x,4d16.8)' ) cby0, cdy0
    write ( *, '(6x,4d16.8)' ) cby1, cdy1
  end subroutine cjy01_test

  subroutine comelp_test(sf)
    type(special_functions), intent(in) :: sf

    real(real64) :: ce, ck, hk
    integer(int32) :: i

    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'COMELP_TEST'
    write ( *, '(a)' ) '  COMELP computes complete elliptic integrals K(K), E(K).'
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '    k         K(k)          E(K)'
    write ( *, '(a)' ) '  ---------------------------------'

    do i = 0, 4
      hk = real(i, kind=real64) / 4.0d0
      call sf%comelp(hk, ck, ce)

      if ( hk /= 1.0d0 ) then
        write ( *, '(2x,f5.2,2f14.6)' ) hk, ck, ce
      else
        write ( *, '(2x,f5.2,3x,a,3x,f14.6)' ) hk, 'Infinity', ce
      end if
    end do
  end subroutine comelp_test

  subroutine hygfx_test(sf)
    type(special_functions), intent(in) :: sf

    real(real64) :: a, b, c, hf, x
    real(real64) :: a_test(4), c_test(4), x_test(3)
    integer(int32) :: i, k, l

    data a_test / -2.5d0, -0.5d0, 0.5d0, 2.5d0 /
    data c_test / -5.5d0, -0.5d0, 0.5d0, 4.5d0 /
    data x_test /  0.25d0, 0.55d0, 0.85d0 /

    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'HYGFX_TEST:'
    write ( *, '(a)' ) '  HYGFX evaluates the hypergeometric function 2F1.'
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '     A              B            C            X             F(A,B,C,X)'

    do l = 1, 3
      x = x_test(l)
      c = 6.7d0
      b = 3.3d0
      write ( *, '(a)' ) ' '
      do i = 1, 4
        a = a_test(i)
        call sf%hygfx(a, b, c, x, hf)
        write ( *, '(4g14.6,g24.16)' ) a, b, c, x, hf
      end do
    end do

    do l = 1, 3
      x = x_test(l)
      write ( *, '(a)' ) ' '
      do k = 1, 4
        c = c_test(k)
        b = 6.7d0
        a = 3.3d0
        call sf%hygfx(a, b, c, x, hf)
        write ( *, '(4g14.6,g24.16)' ) a, b, c, x, hf
      end do
    end do
  end subroutine hygfx_test

  subroutine sphj_test(sf)
    type(special_functions), intent(in) :: sf

    real(real64) :: dj(0:250), sj(0:250), x
    integer(int32) :: k, n, nm, ns

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MSPHJ'
    write ( *, '(a)' ) '  SPHJ evaluates spherical Bessel J functions'

    n = 5
    x = 0.905d0

    if ( n <= 10 ) then
      ns = 1
    else
      ns = 5
    end if

    call sf%sphj(n, x, nm, sj, dj)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   n      x                   jn(x)               jn''(x)'
    write ( *, '(a)' ) ''
    do k = 0, nm, ns
      write ( *, '(1X,I3,3D20.10)' ) k, x, sj(k), dj(k)
    end do

    n = 5
    x = 10.0d0

    if ( n <= 10 ) then
      ns = 1
    else
      ns = 5
    end if

    call sf%sphj(n, x, nm, sj, dj)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   n      x                   jn(x)               jn''(x)'
    write ( *, '(a)' ) ''
    do k = 0, nm, ns
      write ( *, '(1X,I3,3D20.10)' ) k, x, sj(k), dj(k)
    end do
  end subroutine sphj_test

end program TestSpecialFunctions
