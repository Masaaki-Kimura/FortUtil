!=============================================================
!  modfortio_test.f90  (uses ModTestUtil)
!  ModFortIO の動作確認テスト（見やすい出力＋色つき）
!  依存: ModNuclConst.f90, ModFortIO.f90, ModTestUtil.f90
!  ビルド例:
!    gfortran -std=f2008 -O2 ModNuclConst.f90 ModFortIO.f90 ModTestUtil.f90 modfortio_test.f90 -o test_mfio
!=============================================================
program modfortio_test
  use iso_fortran_env, only: real64
  use ModFortIO
  use ModNuclConst
  use ModTestUtil
  implicit none

  type(FileHandler) :: fh
  logical :: ok
  character(:), allocatable :: msg
  character(:), allocatable :: fn_txt, fn_bin
  character(len=256) :: line
  real(real64) :: x, elapsed
  integer :: val_w, val_r, fails

  call t_banner('ModFortIO smoke test with utilities')
  fails = 0

  ! 準備：一時ファイル名
  fn_txt = tmp_filename('mfio_txt', '.txt')
  fn_bin = tmp_filename('mfio_bin', '.bin')
  call t_note('tmp files: '//trim(fn_txt)//', '//trim(fn_bin))

  !-----------------------------------------------------------
  ! [1] Path utilities
  !-----------------------------------------------------------
  call t_section('[1] Path utilities')
  call show_case_c('append_slash', '/a/b', '/a/b/', append_slash('/a/b'))
  call show_case_c('file_path', '/a/b/c.txt', '/a/b/', file_path('/a/b/c.txt'))
  call show_case_c('file_name', '/a/b/c.txt', 'c.txt', file_name('/a/b/c.txt'))
  call show_case_c('file_basename','/a/b/c.tar.gz','c.tar', file_basename('/a/b/c.tar.gz'))
  call show_case_c('file_ext','/a/b/c.tar.gz','.gz', file_ext('/a/b/c.tar.gz'))
  print *, ''

  !-----------------------------------------------------------
  ! [2] file_exists / create_blank_file
  !-----------------------------------------------------------
  call t_section('[2] File existence & blank creation')
  call safe_delete(fn_txt); call safe_delete(fn_bin)
  call assert_bool('file_exists (before)', .false., file_exists(fn_txt), fails)
  ok = fh%create_blank_file(fn_txt)
  call assert_bool('create_blank_file returned', .true., ok, fails)
  call assert_bool('file_exists (after)', .true., file_exists(fn_txt), fails)
  print *, ''

  !-----------------------------------------------------------
  ! [3] Formatted text I/O
  !-----------------------------------------------------------
  call t_section('[3] Formatted text I/O')
  ok = fh%open(fn_txt, status='old', action='readwrite', form='formatted', message=msg)
  if (.not. ok) call t_fatal('open formatted failed: '//trim(msg))

  write(fh%unit,'(A)') 'Line-1: Hello ModFortIO'
  write(fh%unit,'(A)') 'Line-2: test text file'
  ok = fh%close()
  call assert_bool('close(formatted)', .true., ok, fails)

  ok = fh%open(fn_txt, status='old', action='read', form='formatted', message=msg)
  if (.not. ok) call t_fatal('reopen formatted failed: '//trim(msg))
  read(fh%unit,'(A)') line
  call assert_eq_c('read line 1', 'Line-1: Hello ModFortIO', trim(line), fails)
  read(fh%unit,'(A)') line
  call assert_eq_c('read line 2', 'Line-2: test text file', trim(line), fails)
  ok = fh%close()
  call assert_bool('close(readback)', .true., ok, fails)
  print *, ''

  !-----------------------------------------------------------
  ! [4] Binary (stream) I/O
  !-----------------------------------------------------------
  call t_section('[4] Binary (stream) I/O')
  val_w = 123456
  ok = fh%open_binary(fn_bin, status='replace', action='readwrite', message=msg)
  if (.not. ok) call t_fatal('open_binary write failed: '//trim(msg))
  write(fh%unit) val_w
  ok = fh%close()
  call assert_bool('close(binary write)', .true., ok, fails)

  ok = fh%open_binary(fn_bin, status='old', action='read', message=msg)
  if (.not. ok) call t_fatal('open_binary read failed: '//trim(msg))
  read(fh%unit) val_r
  ok = fh%close()
  call assert_bool('close(binary read)', .true., ok, fails)
  call assert_eq_i('binary value', val_w, val_r, fails)
  print *, ''

  !-----------------------------------------------------------
  ! [5] IEEE NaN
  !-----------------------------------------------------------
  call t_section('[5] IEEE NaN')
  x = ieee_nan()
  call assert_true('isnan(ieee_nan())', (x /= x), fails)
  print *, ''

  !-----------------------------------------------------------
  ! [6] text_color demo
  !-----------------------------------------------------------
  call t_section('[6] text_color demo')
  print *, text_color('g', '  Success: green text sample')
  print *, text_color('r', '  Error:   red text sample')
  print *, text_color('b', '  Info:    blue text sample')
  call t_note('色が出ない環境では modfortio_enable_text_color=.false. を推奨')
  print *, ''

  !-----------------------------------------------------------
  ! [7] Timer / RNG / argv (from ModTestUtil)
  !-----------------------------------------------------------
  call t_section('[7] Timer / RNG / argv (from ModTestUtil)')
  call tic()
  call seed_rng()
  call dump_argv()
  elapsed = toc()   ! 関数なので call ではなく代入で受ける
  call t_info('elapsed(s) = '//trim(adjustl(to_str_real(elapsed))))
  print *, ''

  !-----------------------------------------------------------
  ! [8] Cleanup
  !-----------------------------------------------------------
  call t_section('[8] Cleanup')
  call safe_delete(fn_txt); call safe_delete(fn_bin)
  call assert_bool('file_exists(txt) after delete', .false., file_exists(fn_txt), fails)
  call assert_bool('file_exists(bin) after delete', .false., file_exists(fn_bin), fails)
  print *, ''

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

contains
  ! 簡易：real を文字列化（I/O依存なので pure にはしない）
  function to_str_real(x) result(s)
    real(real64), intent(in) :: x
    character(len=64) :: s
    write(s,'(ES12.5)') x
  end function to_str_real
end program modfortio_test
