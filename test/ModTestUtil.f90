!=============================================================
!  ModTestUtil.f90  (fixed)
!  汎用テストユーティリティ（見やすい出力、assert群、I/O補助、計時 等）
!  依存: iso_fortran_env, ModFortIO（色付け）
!  ※ 先に ModFortIO をコンパイルしてください。
!=============================================================
module ModTestUtil
  use iso_fortran_env, only: int64, real64
  use ModFortIO,      only: text_color, file_exists
  implicit none
  private
  public :: t_banner, t_section, t_info, t_note, t_fatal
  public :: t_okmark, t_ngmark
  public :: assert_eq_i, assert_eq_r, assert_eq_c, assert_true, assert_false, assert_bool
  public :: show_case_c, show_bool
  public :: tmp_filename, safe_delete, write_file_lines, read_file_lines
  public :: tic, toc
  public :: seed_rng, dump_argv

  integer(int64), save :: tic_count_ = 0_int64   ! ← 先頭アンダースコア禁止なのでリネーム
contains
  !------------- 表示系 -------------
  subroutine t_banner(title)
    character(*), intent(in) :: title
    print *, repeat('=', 64)
    print *, text_color('c', '  '//trim(title))
    print *, repeat('=', 64)
    print *, ''
  end subroutine t_banner

  subroutine t_section(title)
    character(*), intent(in) :: title
    print *, text_color('b', trim(title))
  end subroutine t_section

  subroutine t_info(msg)
    character(*), intent(in) :: msg
    print *, text_color('c','  [INFO] ')//trim(msg)
  end subroutine t_info

  subroutine t_note(msg)
    character(*), intent(in) :: msg
    print *, text_color('y','  [NOTE] ')//trim(msg)
  end subroutine t_note

  subroutine t_fatal(msg)
    character(*), intent(in) :: msg
    print *, text_color('r','  [FATAL] ')//trim(msg)
    error stop 'fatal'
  end subroutine t_fatal

  subroutine t_okmark()
    print *, '   ', text_color('g','OK')
  end subroutine t_okmark

  subroutine t_ngmark()
    print *, '   ', text_color('r','NG')
  end subroutine t_ngmark

  !------------- アサート -------------
  subroutine assert_eq_i(label, expected, got, fails)
    character(*), intent(in) :: label
    integer,      intent(in) :: expected, got
    integer,      intent(inout) :: fails
    write(*,'(A,1X,A,1X,A,1X,I0,1X,A,1X,I0)') '  ', trim(label), &
         'expected=', expected, 'got=', got
    if (got == expected) then
      call t_okmark()
    else
      call t_ngmark(); fails = fails + 1
    end if
  end subroutine assert_eq_i

  subroutine assert_eq_r(label, expected, got, tol, fails)
    character(*), intent(in) :: label
    real(real64), intent(in) :: expected, got, tol
    integer,      intent(inout) :: fails
    write(*,'(A,1X,A,1X,A,1X,ES12.5,1X,A,1X,ES12.5,1X,A,1X,ES10.3)') &
         '  ', trim(label), 'expected=', expected, 'got=', got, 'tol=', tol
    if (abs(got-expected) <= tol) then
      call t_okmark()
    else
      call t_ngmark(); fails = fails + 1
    end if
  end subroutine assert_eq_r

  subroutine assert_eq_c(label, expected, got, fails)
    character(*), intent(in) :: label, expected, got
    integer,      intent(inout) :: fails
    write(*,'(A,1X,A,1X,A,1X,A,1X,A)') '  ', trim(label), &
         'expected=', '"'//trim(expected)//'"', 'got=', '"'//trim(got)//'"'
    if (trim(got) == trim(expected)) then
      call t_okmark()
    else
      call t_ngmark(); fails = fails + 1
    end if
  end subroutine assert_eq_c

  subroutine assert_true(label, expr, fails)
    character(*), intent(in) :: label
    logical,      intent(in) :: expr
    integer,      intent(inout) :: fails
    write(*,'(A,1X,A,1X,A,1X,A)') '  ', trim(label), 'expected=','T', 'got= '//merge('T','F',expr)
    if (expr) then
      call t_okmark()
    else
      call t_ngmark(); fails = fails + 1
    end if
  end subroutine assert_true

  subroutine assert_false(label, expr, fails)
    character(*), intent(in) :: label
    logical,      intent(in) :: expr
    integer,      intent(inout) :: fails
    write(*,'(A,1X,A,1X,A,1X,A)') '  ', trim(label), 'expected=','F', 'got= '//merge('T','F',expr)
    if (.not. expr) then
      call t_okmark()
    else
      call t_ngmark(); fails = fails + 1
    end if
  end subroutine assert_false

  subroutine assert_bool(label, expected, got, fails)
    character(*), intent(in) :: label
    logical,      intent(in) :: expected, got
    integer,      intent(inout) :: fails
    write(*,'(A,1X,A,1X,A,1X,A)') '  ', trim(label), 'expected=', merge('T','F',expected), &
                                   'got= '//merge('T','F',got)
    if (got .eqv. expected) then
      call t_okmark()
    else
      call t_ngmark(); fails = fails + 1
    end if
  end subroutine assert_bool

  !------------- 表示ユーティリティ -------------
  subroutine show_case_c(label, input, expected, got)
    character(*), intent(in) :: label, input, expected, got
    print *, '  '//trim(label)//'  input:    "'//trim(input)//'"'
    print *, '                     expected: "'//trim(expected)//'"'
    print *, '                     got:      "'//trim(got)//'"'
    if (trim(got) == trim(expected)) then
      call t_okmark()
    else
      call t_ngmark()
    end if
  end subroutine show_case_c

  subroutine show_bool(label, expected, got)
    character(*), intent(in) :: label
    logical,      intent(in) :: expected, got
    print *, '  '//trim(label)//'  expected:', merge('T','F',expected), '  got:', merge('T','F',got)
    if (got .eqv. expected) then
      call t_okmark()
    else
      call t_ngmark()
    end if
  end subroutine show_bool

  !------------- I/O ユーティリティ -------------
  function tmp_filename(prefix, ext) result(fn)
    ! PUREは外す（system_clock使用）
    character(*), intent(in) :: prefix
    character(*), intent(in), optional :: ext
    character(:), allocatable :: fn
    integer :: clk, pid
    character(64) :: buf
    call system_clock(count=clk)
    pid = abs(1103515245*clk + 12345)   ! 緩い一意化
    write(buf,'(A,"_",I0)') trim(prefix), pid
    if (present(ext)) then
      fn = trim(buf)//trim(ext)
    else
      fn = trim(buf)//'.tmp'
    end if
  end function tmp_filename

  subroutine safe_delete(fn)
    character(*), intent(in) :: fn
    if (file_exists(fn)) then
      open(unit=99, file=fn, status='old')
      close(99, status='delete')
    end if
  end subroutine safe_delete

  subroutine write_file_lines(fn, lines, ios, msg)
    character(*), intent(in)  :: fn
    character(*), intent(in)  :: lines(:)
    integer,      intent(out) :: ios
    character(:), allocatable, intent(out) :: msg
    integer :: u, i
    ios = 0; msg=''
    open(newunit=u, file=fn, status='replace', action='write', iostat=ios)
    if (ios/=0) then
      msg='write_file_lines: open failed'; return
    end if
    do i=1,size(lines)
      write(u,'(A)', iostat=ios) trim(lines(i))
      if (ios/=0) then
        msg='write_file_lines: write failed'; exit
      end if
    end do
    close(u, iostat=ios)
    if (ios/=0) msg='write_file_lines: close failed'
  end subroutine write_file_lines

  subroutine read_file_lines(fn, lines, ios, msg, maxn)
    character(*), intent(in)  :: fn
    character(:), allocatable, intent(out) :: lines(:)
    integer,      intent(out) :: ios
    character(:), allocatable, intent(out) :: msg
    integer,      intent(in),  optional :: maxn
    integer :: u, nmax, n, stat
    character(len=1024) :: buf
    integer :: L
    ios = 0; msg=''; n=0
    if (present(maxn)) then
      nmax = maxn
    else
      nmax = 100000
    end if
    ! 要素長を明示して0配列を確保（deferred-length対応）
    allocate(character(len=1) :: lines(0))
    open(newunit=u, file=fn, status='old', action='read', iostat=ios)
    if (ios/=0) then
      msg='read_file_lines: open failed'; return
    end if
    do
      read(u,'(A)', iostat=stat) buf
      if (stat /= 0) exit
      if (n >= nmax) exit
      n = n + 1
      L = len_trim(buf)
      call push_line(lines, buf(1:L))
    end do
    close(u)

  contains
  subroutine push_line(a, s)
    character(:), allocatable, intent(inout) :: a(:)
    character(*), intent(in) :: s
    character(:), allocatable :: tmp(:)
    integer :: m, newlen
    m = size(a)
    newlen = max(len(a), len_trim(s))   ! len(a): 要素長（aは事前にallocate済み）
    allocate(character(len=newlen) :: tmp(m+1))
    tmp = ''                             ! ★ これを追加（未初期化警告の解消）
    if (m>0) tmp(1:m) = a
    tmp(m+1) = s
    call move_alloc(tmp, a)
  end subroutine push_line  
  end subroutine read_file_lines


  !------------- 計時 -------------
  subroutine tic()
    integer(int64) :: c
    call system_clock(count=c)
    tic_count_ = c
  end subroutine tic

  real(real64) function toc() result(elapsed_sec)
    integer(int64) :: c, rate
    call system_clock(count=c, count_rate=rate)
    if (tic_count_ <= 0_int64) then
      elapsed_sec = 0.0_real64
    else
      elapsed_sec = real(c - tic_count_, kind=real64) / real(rate, kind=real64)
    end if
    print *, text_color('m', '  [TIMER] elapsed (s) = '), elapsed_sec
  end function toc

  !------------- RNG -------------
  subroutine seed_rng(seed)
    integer, intent(in), optional :: seed
    integer, allocatable :: s(:)
    integer :: n, i, base
    call random_seed(size=n)
    allocate(s(n))
    if (present(seed)) then
      base = seed
    else
      call system_clock(count=base)
    end if
    do i=1,n
      s(i) = ieor(base, 7919*i)  ! ← ieor に修正（ixorは非標準）
    end do
    call random_seed(put=s)
  end subroutine seed_rng


  subroutine dump_argv()
    integer :: n, i
    character(len=1024) :: buf
    call get_command_argument(0, buf)
    print *,'[ARGV] program: ', trim(buf)
    n = command_argument_count()      ! ★ サブルーチン呼び → 関数呼びに修正
    do i=1,n
      call get_command_argument(i, buf)
      print *,'[ARGV] arg', i, ': ', trim(buf)
    end do
  end subroutine dump_argv

end module ModTestUtil
