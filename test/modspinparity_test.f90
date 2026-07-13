!=============================================================
!  modspinparity_test.f90  (uses ModTestUtil)
!=============================================================
program modspinparity_test
  use ModFortIO
  use ModNuclConst
  use ModSpinParity
  use ModTestUtil
  implicit none

  integer :: fails
  type(SpinParity) :: a, b, c, s_only
  logical :: ok
  character(:), allocatable :: s

  call t_banner('ModSpinParity pretty test')
  fails = 0

  ! [1] Constructors & sp_str()
  call t_section('[1] Constructors & sp_str')
  a = SpinParity('2+')
  call assert_eq_i('ctor "2+" -> j', 2, a%j, fails)
  call assert_eq_i('ctor "2+" -> p', POS, a%p, fails)
  s = a%str()
  call assert_eq_c('sp_str("2+")', '2+', s, fails)

  b = SpinParity(2, NEG)
  call assert_eq_i('ctor (2,NEG) -> j', 2, b%j, fails)
  call assert_eq_i('ctor (2,NEG) -> p', NEG, b%p, fails)
  s = b%str()
  call assert_eq_c('sp_str(2,NEG)', '2-', s, fails)

  c = SpinParity(0, POS)
  call assert_eq_i('ctor (0,POS) -> j', 0,  c%j, fails)
  call assert_eq_i('ctor (0,POS) -> p', POS,c%p, fails)
  call assert_eq_i('dim = j+1', c%j+1, c%dim, fails)
  print *, ''

  ! [2] m array
  call t_section('[2] m array mapping')
  call assert_eq_i('a%dim (j=2)', 3, a%dim, fails)
  call assert_eq_i('a%m(1) = -2', -2, a%m(1), fails)
  call assert_eq_i('a%m(2) =   0',  0, a%m(2), fails)
  call assert_eq_i('a%m(3) =  +2', +2, a%m(3), fails)
  print *, ''

  ! [3] is_natural()  (π = (-1)^J)
  call t_section('[3] is_natural()')
  a = SpinParity(0, POS); call assert_true ('nat(0+)', a%is_natural(), fails)
  a = SpinParity(2, NEG); call assert_true ('nat(2-)', a%is_natural(), fails)
  a = SpinParity(2, POS); call assert_false('nat(2+)', a%is_natural(), fails)
  a = SpinParity(4, POS); call assert_true ('nat(4+)', a%is_natural(), fails)
  a = SpinParity(4, NEG); call assert_false('nat(4-)', a%is_natural(), fails)
  print *, ''

  ! [4] is_triangle()
  call t_section('[4] is_triangle()')
  a = SpinParity(2, POS);  b = SpinParity(2, POS);  c = SpinParity(0, POS)
  ok = a%is_triangle(b, c); call assert_true ('triangle: (2+ ⊗ 2+ -> 0+)', ok, fails)

  b = SpinParity(2, NEG)
  ok = a%is_triangle(b, c); call assert_false('triangle: (2+ ⊗ 2- -> 0+)', ok, fails)

  a = SpinParity(2, NEG);  b = SpinParity(2, NEG);  c = SpinParity(4, POS)
  ok = a%is_triangle(b, c); call assert_true ('triangle: (1- ⊗ 1- -> 2+)', ok, fails)
  print *, ''

  ! [5] spin-only（ダミーで netural parity を付与）
  call t_section('[5] spin-only with dummy natural parity')
  s_only = SpinParity(6, NEG) ! J=3, π=-
  call assert_eq_i('spin-only j (2J)', 6, s_only%j, fails)
  call assert_true('spin-only natural?', s_only%is_natural(), fails)
  a = SpinParity(2, POS); b = SpinParity(4, NEG); c = s_only
  ok = a%is_triangle(b, c)
  call assert_true('triangle: (1+ ⊗ 2+ -> 3-)', ok, fails)
  print *, ''

  ! 結果
  if (fails == 0) then
    print *, repeat('=', 64)
    print *, text_color('g','ALL TESTS PASSED')
    print *, repeat('=', 64)
  else
    print *, repeat('-', 64)
    write(*,'(A,I0)') text_color('r','FAILED TESTS: '), fails
    print *, repeat('-', 64)
    error stop 'some tests failed'
  end if
end program modspinparity_test
