!> @file ModFortIO.f90
!! @brief Module of I/O utilities in Fortran
!!
module ModFortIO
  use iso_fortran_env, only: real64, input_unit, output_unit, error_unit
  implicit none
  private

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! public module variables
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! control text color
  logical, save, public :: mfio_enable_text_color = .true. ! enable text color 
  character(len=1), public :: mfio_error_color = 'r' ! set error text color

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! public constants
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! default size of character buffer
  integer, parameter, public :: MAXBUF = 512
  ! reserved I/O unit numbers
  integer, parameter, public :: STDERR = error_unit
  integer, parameter, public :: STDIN  = input_unit
  integer, parameter, public :: STDOUT = output_unit
  ! ASCII for TAB, CR, LF
  character(len=1), parameter, public :: ASCII_TAB = char(9)
  character(len=1), parameter, public :: ASCII_LF  = char(10)
  character(len=1), parameter, public :: ASCII_CR  = char(13)
  character(len=1), parameter, public :: ASCII_ESC = char(27)


  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! public class
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public :: FileHandler    !< class for file I/O

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !! public functions/subroutines
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! change text color in terminal
  public :: text_color     !< set text color

  ! file name utility
  public :: append_slash   !< append '/' to the input; append_slash('/foo/bar') -> '/foo/bar/'
  public :: file_path      !< get path; file_path('/foo/bar/sample.tar.gz') -> '/foo/bar/'
  public :: file_name      !< get file name; file_name('/foo/bar/sample.tar.gz') -> 'sample.tar.gz'
  public :: file_basename  !< get basename; file_basename('/foo/bar/sample.tar.gz') -> 'sample.tar'
  public :: file_ext       !< get extension; file_ext('/foo/bar/sample.tar.gz') -> '.gz'
  public :: file_exists    !< return .true. if a file exists; file_exists('/foo/bar/sample.tar.gz')

  ! return NaN_real64
  public :: ieee_nan       !< return a NaN (not-a-number) with real(real64) type

  ! command line argument parser
  public :: get_arg        !! get command line argument for given position
  public :: next_arg       !! get next command line argument

  ! conversion to string (int -> string, real -> string, logical -> string)
  public :: to_string      !! convert various types to string
  public :: to_int         !! convert string to integer
  public :: to_real        !! convert string to real

  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !>
  !> @brief Class for file I/O\
  !! integer::unit : unit number associated to file fh%name\
  !! character(:), allocatable :: name : file name associated to unit fh%unit\
  !! logical open(fn,status,action,position,form,message): open a file and set fh%unit\
  !! logical open_binary(fn,status,action,position,message): open a binary file and set fh%unit\
  !! logical close(): close the file associated to fh%unit and set fh%unit to error_unit\
  !! logical create_blank_file(fn): create a blank file in a new mode\
  !! final filehandler_delete: destructor
  type FileHandler
    integer :: unit = error_unit  ! initialize with error_unit
    character(:), allocatable :: name

    contains
    procedure, public :: open  => filehandler_open
    procedure, public :: open_binary => filehandler_open_binary
    procedure, public :: close => filehandler_close
    procedure, nopass, public :: create_blank_file => filehandler_create_blank_file

    final :: delete_filehandler
  end type FileHandler

  !>
  !> @brief constructor for FileHandler
  interface FileHandler
    module procedure init_filehandler
  end interface FileHandler


  !> @brief converter to string
  interface to_string
    module procedure to_string_int, to_string_real, to_string_logical
  end interface to_string


  contains

  !> @brief constructor for FileHandler
  !! @return FileHandler object
  pure type(FileHandler) function init_filehandler()
    implicit none
    init_filehandler%unit = error_unit
  end function init_filehandler

  !> @brief destructor for FileHandler
  !! @param this: FileHandler object to be deleted
  subroutine delete_filehandler(this)
    implicit none
    ! argument
    type(FileHandler), intent(inout) :: this ! myself
    ! local variables
    logical :: dammy
    
    ! body
    if(this%unit /= error_unit) dammy = this%close() ! close file if it is opened

  end subroutine delete_filehandler
    

  !> @brief Open a file with a new unit number and set fh%unit to the unit number\
  !! @notes If the file is unformatted, access is forced to 'stream'.\
  !! @notes If the file is formatted, access is forced to 'sequential'.\
  !! @notes If the file is already opened, return false and show an error message.\
  !! @return .true. on success
  !! @param fn: file name
  !! @param status: file open status ('old','new','replace','scratch','unknown')
  !!                default: 'unknown'
  !! @param action: file open action ('read','write','readwrite')
  !!                default: 'readwrite'
  !! @param form: file form ('formatted','unformatted')
  !!              default: 'formatted'
  !! @param position: file position ('rewind','append','asis')
  !!                  default: 'rewind'
  !! @param message: (optional) error message on failure
  logical function filehandler_open(this,fn,status,action,position,form,message)
    implicit none
    ! argument
    class(FileHandler), intent(inout) :: this
    character(len=*), intent(in) :: fn
    character(len=*), intent(in), optional :: status
    character(len=*), intent(in), optional :: action
    character(len=*), intent(in), optional :: form
    character(len=*), intent(in), optional :: position
    character(:), allocatable, intent(out), optional :: message
    ! local variables 
    integer :: fstat
    logical :: is_opened
    character(:), allocatable :: stat, act, pos, frm
    character(len=MAXBUF) :: iomsg
    ! body

    filehandler_open = .false. ! default rerun value

    if(this%unit /= error_unit) then    ! check if this%unit is already opened
      write(error_unit,'(A)') &
      text_color(mfio_error_color,'filehandler_open: this instance is already associated to an opened file: '//trim(this%name))
      return
    end if

    ! Check if the file 'fn' is already opened
    inquire(file=fn, opened=is_opened)
    if(is_opened) then    ! if the file is already opened, set unit number and return false.
      write(error_unit,'(A)') &
      text_color(mfio_error_color,'filehandler_open:')//' '//trim(fn)//' is already opened.'
      return
    end if

    ! set optional arguments
    stat='unknown'   ; if(present(status)) stat = status
    act='readwrite'  ; if(present(action)) act = action 
    pos='rewind'     ; if(present(position)) pos = position 
    frm='formatted'  ; if(present(form)) frm = form 


    ! Open file with new unit number
    if(trim(frm) == 'unformatted') then      ! for unformatted file, force access to 'stream'
      open(newunit=this%unit,file=fn,status=stat,form=frm,access='stream',action=act,position=pos,iostat=fstat,iomsg=iomsg)
    else ! for formatted file
      open(newunit=this%unit,file=fn,status=stat,form=frm,access='sequential',action=act,position=pos,iostat=fstat,iomsg=iomsg)
    end if

    ! check open status and return false on error
    if(fstat/=0) then
      if(present(message)) message = trim(adjustl(iomsg))
      write(error_unit,'(A)') text_color(mfio_error_color,'filehandler_open: ')//' cannot open '//trim(fn)
      this%unit = error_unit
      return
    end if
    ! set file name
    this%name = trim(adjustl(fn))
    filehandler_open = .true.
  end function filehandler_open


  !> @brief Open a binary file with a new unit number and set fh%unit to the unit number\
  !! @notes A wrapper function of filehandler_open with form='unformatted'.\
  !! @notes If the file is already opened, return false and show an error message.
  !! @return .true. on success
  !! @param fn: file name
  !! @param status: file open status ('old','new','replace','scratch','unknown')
  !!                default: 'unknown'
  !! @param action: file open action ('read','write','readwrite')
  !!                default: 'readwrite'
  !! @param position: file position ('rewind','append','asis')
  !!                  default: 'rewind'
  !! @param message: (optional) error message on failure
  logical function filehandler_open_binary(this,fn,status,action,position,message)
    implicit none
    ! argument
    class(FileHandler), intent(inout) :: this
    character(len=*), intent(in) :: fn
    character(len=*), intent(in), optional :: status
    character(len=*), intent(in), optional :: action
    character(len=*), intent(in), optional :: position
    character(:), allocatable, intent(out), optional :: message
    ! local variables 
    logical :: is_opened
    character(:), allocatable :: stat, act, pos, iomsg
    ! body

    filehandler_open_binary = .false. ! default rerun value
    stat='unknown' ; if(present(status)) stat = status    ! set status
    act='readwrite'; if(present(action)) act = action     ! set action
    pos='rewind'   ; if(present(position)) pos = position ! set position

    ! Check if the file 'fn' is already opened
    inquire(file=fn, opened=is_opened)
    if(is_opened) then    ! if the file is already opened, set unit number and return false.
      write(error_unit,'(A)') &
      text_color(mfio_error_color,'filehandler_open_binary: ')//trim(fn)//' is already opened.'
      return
    end if

    ! Open file with new unit number
    filehandler_open_binary = this%open(fn,status=stat,action=act,position=pos, &
                                           form='unformatted',message=iomsg)
    if(.not. filehandler_open_binary) then
      write(error_unit,'(A)') &
      text_color(mfio_error_color,'filehandler_open_binary: ')//' cannot open '//trim(fn)
      if(present(message)) message = trim(adjustl(iomsg))
      return
    end if
  end function filehandler_open_binary

  !> @brief close the file associated to fh%unit and set fh%unit to error_unit
  !! @return .true. on success
  logical function filehandler_close(this)
    implicit none
    ! argument
    class(FileHandler), intent(inout) :: this
    ! local variables
    integer :: ios
    logical :: is_file_opened
    ! body
    filehandler_close = .false.    ! set default return value

    ! check if the file with name of this%name is already opened
    inquire(unit=this%unit, opened=is_file_opened)
    if(.not. is_file_opened) then ! if the file is not opened, show a warning and return false
      write(error_unit,'(A)') text_color(mfio_error_color,'filehandler_close:')//' file is not opened.'
      return
    end if

    close(this%unit,iostat=ios)    ! close the file
    if(ios /= 0) then               ! on error, show message and return
      write(error_unit,'(A)') text_color(mfio_error_color,'filehandler_close:')//' cannot close the file: '//trim(this%name)
      return
    end if
    this%unit = error_unit  ! set to error_unit
    filehandler_close = .true. ! return true on success
  end function filehandler_close
 

  !> @brief create a blank file in a new mode.
  !! This does not set this%name
  !! @return .true. on success, .false. if the file already exists
  !! @param fn: file name
  logical function filehandler_create_blank_file(fn)
    implicit none
    ! arguments
    character(len=*), intent(in) :: fn
    ! local variables
    integer :: unit,ios
    
    ! body
    filehandler_create_blank_file = .false. ! default return value
    open(newunit=unit,file=trim(fn),status='new',action='write',iostat=ios) ! create a blank file
    if (ios /= 0) then ! On error, show message and return
      write(error_unit,"(A)") text_color(mfio_error_color,"filehandler_create_blank_file:")//" cannot create a new file: "//trim(fn)
      return
    end if

    close(unit, iostat=ios) ! close the file
    if (ios /= 0) then ! On error, show message and return
      write(error_unit,"(A)") text_color(mfio_error_color,"filehandler_create_blank_file:")//" cannot close the file: "//trim(fn)
      return
    end if

    filehandler_create_blank_file = .true.  ! return true on success
    return

  end function filehandler_create_blank_file


  !> @brief set text color
  !! @return colored text
  !! @param color: color code ('r':red, 'g':green, 'y':yellow, 'b':blue, 'm':magenta, 'c':cyan)
  !! @param src: source text
  pure function text_color(color,src)
    implicit none
    ! return value
    character(:), allocatable :: text_color
    ! arguments
    character(len=*), intent(in) :: color
    character(len=*), intent(in) :: src
    ! local variables
    character(len=2) :: ColorCode

    ! body
    if(.not. mfio_enable_text_color) then ! if text color is disabled, return src
      text_color = trim(adjustl(src))
      return
    end if
    ! set color code 
    select case (color)
    case('r') ! red
      ColorCode = '31'
    case('g') ! green
      ColorCode = '32'
    case('y') ! yellow
      ColorCode = '33'
    case('b') ! blue
      ColorCode = '34'
    case('m') ! magenta
      ColorCode = '35'
    case('c') ! cyan
      ColorCode = '36'
    case default ! default
      ColorCode = '0'
    end select

    text_color = ASCII_ESC//'['//ColorCode//'m'//trim(adjustl(src))//ASCII_ESC//'[0m'
  
  end function text_color


  !////////////////////////////////////////////////////////////////////////////////////////////////
  ! functions for file name formatting
  !////////////////////////////////////////////////////////////////////////////////////////////////
  
  !> @brief append '/' to the input
  !! @return string with '/' at the end
  !! @param src: input string
  pure function append_slash(src)
    implicit none
    ! return value
    character(:), allocatable :: append_slash
    ! arguments
    character(len=*),intent(in) :: src

    ! body
    append_slash = trim(adjustl(src))
    ! add slash if the last character is not '/'
    if(len(append_slash) /= index(append_slash,'/',back=.true.)) append_slash = append_slash//'/'

  end function append_slash
  

  !> @brief Get path to file including the last slash, \
  !! e.g. file_path('/foo/bar/sample.tar.gz') -> /foo/bar/
  !! @return path to file including the last slash
  !! @param src: input string
  pure function file_path(src)
    implicit none
    ! return value
    character(:), allocatable :: file_path
    ! arguments
    character(len=*), intent(in) :: src
    ! local variables
    integer :: slash
    
    ! body
    ! copy src
    file_path = trim(adjustl(src))
    ! get position of the last slash
    slash = index(file_path,'/',back=.true.)
    
    ! return path
    if(slash <= 0) then ! if no slash is found, return current position
      file_path = './' ! current position
    else ! if slash is found, return text before the last slash 
      file_path = file_path(:slash)
    end if

  end function file_path  


  !> @brief Get file name, \
  !! e.g. file_name('/foo/bar/sample.tar.gz') -> sample.tar.gz
  !! @return file name
  !! @param src: input string
  pure function file_name(src)
    implicit none
    ! return value
    character(:), allocatable :: file_name
    ! arguments
    character(len=*), intent(in) :: src
    ! local variables
    integer :: slash

    ! body
    file_name = trim(adjustl(src))
    ! get position of the last slash
    slash = index(file_name,'/',back=.true.)

    ! retrun src after the last slash
    file_name = file_name(slash+1:)

  end function file_name


  !> @brief Get basename of a file,\
  !! e.g. file_basename('/foo/bar/sample.tar.gz') -> sample.tar
  !! @return basename
  !! @param src: input string
  pure function file_basename(src)
    implicit none
    ! return value
    character(:), allocatable :: file_basename
    ! arguments
    character(len=*), intent(in) :: src
    ! local variables
    integer :: slash, dot

    ! body
    file_basename = trim(adjustl(src))
    ! get position of the last slash and the last dot
    slash = index(file_basename,'/',back=.true.)
    dot = index(file_basename,'.',back=.true.)
    ! get basename
    if(dot <= slash) then ! if no dot after slash (or no-slash no-dot), return src after slash
      file_basename = file_basename(slash+1:)
    else
      file_basename = file_basename(slash+1:dot-1) ! return src between slash and dot
    end if

  end function file_basename
  
  
  !> @brief Get extension of a file,\
  !! e.g. file_ext('/foo/bar/sample.tar.gz') -> .gz
  !! @return extension
  !! @param src: input string
  pure function file_ext(src)
    implicit none
    ! return value
    character(:), allocatable :: file_ext
    ! arguments
    character(len=*), intent(in) :: src
    ! local variables
    integer :: slash, dot

    ! body
    file_ext = trim(adjustl(src))
    ! get position of the last slash and the last dot
    slash = index(file_ext,'/',back=.true.)
    dot = index(file_ext,'.',back=.true.)
    
    ! get extension
    if(dot <= 0 .or. dot <= slash) then ! no period after slash, return empty string
      file_ext = ''
    else
      file_ext = trim(src(index(file_ext,'.',back=.true.):)) ! return src after the last dot
    end if
  end function file_ext
  

  !> @brief Inquire if the file 'fn' exists
  !> @return .true. if the file 'fn' exists, .false. otherwise
  !> @param fn: file name
  logical function file_exists(fn)
    implicit none
    ! arguments
    character(len=*), intent(in) :: fn
    ! local variables

    ! body
    inquire(file=fn,exist=file_exists) ! inquire if the file 'fn' exists

  end function file_exists


  !> @brief return a NaN (not-a-number) with real(real64) type
  !! @return NaN with real(real64) type
  pure real(real64) function ieee_nan()
    use, intrinsic :: ieee_arithmetic
    implicit none
    ! body
    ieee_nan = ieee_value(0.0_real64, ieee_quiet_nan)
  end function ieee_nan


  !> @brief Get command line argument for given position
  !! @return command line argument as a string
  !!         on error, return an empty string
  !! @param  n position of argument to get
  function get_arg(n)
    implicit none
    ! return value
    character(:), allocatable :: get_arg
    ! arguments
    integer, intent(in) :: n
    ! local variables
    integer :: ios
    character(len=MAXBUF) :: buf

    get_arg = ''
    call get_command_argument(n, buf, status=ios)
    if(ios == 0) get_arg = trim(adjustl(buf)) 
  end function  


  !> @brief Get next command line argument
  !! @return next command line argument as a string.\
  !!         On error, return an empty string
  function next_arg()
    implicit none
    ! return value
    character(:), allocatable :: next_arg
    ! local variables
    integer :: ios
    integer, save :: p = 1
    character(len=MAXBUF) :: buf

    next_arg = ''
    call get_command_argument(p, buf, status=ios)
    if(ios /= 0) return ! return empty string on error
    next_arg = trim(adjustl(buf)) 
    p = p + 1

  end function  


  !> @brief convert integer to string
  !! @return string converted from integer
  !! @param val: integer value to be converted
  !! @param fmt: (optional) format for conversion
  pure function to_string_int(val,fmt)
    implicit none
    ! return value
    character(:), allocatable :: to_string_int
    ! arguments
    integer, intent(in) :: val
    character(len=*), intent(in), optional :: fmt ! format
    ! local variables
    integer :: ios
    character(len=MAXBUF) :: buf
    ! body 
    write(buf,'(I0)',iostat=ios) val ! convert to string with default format
    if(present(fmt)) write(buf,trim(fmt),iostat=ios) val ! convert to string with specified format
    if(ios /= 0) then
      to_string_int = text_color(mfio_error_color,'to_string_int:')//' conversion error.'
      return
    end if
    to_string_int = trim(adjustl(buf))

  end function to_string_int


  !> @brief convert real(real64) to string
  !! @return string converted from real(real64  )
  !! @param val: real(real64) value to be converted
  !! @param fmt: (optional) format for conversion
  pure function to_string_real(val,fmt)
    implicit none
    ! return value
    character(:), allocatable :: to_string_real
    ! arguments
    real(real64), intent(in) :: val
    character(len=*), intent(in), optional :: fmt ! format
    ! local variables
    integer :: ios
    character(len=MAXBUF) :: buf
    ! body 
    write(buf,'(G0)',iostat=ios) val ! convert to string with default format
    if(present(fmt)) write(buf,trim(fmt),iostat=ios) val ! convert to string with specified format
    if(ios /= 0) then
      to_string_real = text_color(mfio_error_color,'to_string_real:')//' conversion error.'
      return
    end if
    to_string_real = trim(adjustl(buf)) 

  end function to_string_real


  !> @brief convert logical to string
  !! @return string converted from logical
  !! @param val: logical value to be converted
  !! @param lang: (optional) language of the output string ('Fortran','C','Python')
  pure function to_string_logical(val,lang)
    implicit none
    ! return value
    character(:), allocatable :: to_string_logical
    ! arguments
    logical, intent(in) :: val
    character(len=*), intent(in), optional :: lang ! language (default: fortran)
    ! body
    
    ! default case (fortran)
    to_string_logical = '.FALSE.'
    if(val) to_string_logical = '.TRUE.'
    
    ! lang is given
    if(present(lang)) then
      select case (trim(adjustl(lang)))
      case('Fortran','fortran','FORTRAN') ! Fortran
        to_string_logical = '.FALSE.'
        if(val) to_string_logical = '.TRUE.'
      case('C','c') ! C language
        to_string_logical = 'false'
        if(val) to_string_logical = 'true'
      case('Python','python','PYTHON') ! Python
        to_string_logical = 'False'
        if(val) to_string_logical = 'True'
      case default
        to_string_logical = &
        text_color(mfio_error_color,'to_string_logical:')//' unknown language: '//trim(adjustl(lang))
      end select
    end if

  end function to_string_logical

  !> @brief convert string to integer
  !! @return integer converted from string
  pure integer function to_int(str)
    implicit none
    ! arguments
    character(len=*), intent(in) :: str
    ! local variables
    character(:), allocatable :: buf
    ! body
    buf = trim(adjustl(str)) ! copy to internal buffer to be a pure function
    read(buf,*) to_int       ! and read from internal buffer

  end function to_int


  !> @brief convert string to real(real64)
  !! @return real(real64) converted from string
  pure real(real64) function to_real(str)
    implicit none
    ! arguments
    character(len=*), intent(in) :: str
    ! local variables
    character(:), allocatable :: buf
    ! body
    buf = trim(adjustl(str)) ! copy to internal buffer to be a pure function
    read(buf,*) to_real      ! and read from internal buffer

  end function to_real


end module ModFortIO

  