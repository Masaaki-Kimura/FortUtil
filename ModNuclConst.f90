!> @file ModNuclConst.f90
!! @brief Module for nuclear constants
module ModNuclConst
  use iso_fortran_env, only: real64
  implicit none
  public

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Definition of mathematical and physical constants
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !> Imaginary unit
  complex(real64), parameter :: I_IMAG=cmplx(0.0_real64,1.0_real64,kind=real64)

  !> Value of PI
  real(real64), parameter :: PI = acos(-1.0_real64)

  !> Fine structure constant
  real(real64), parameter :: FS_CONST = 0.0072973525643_real64

  !> Proton mass (MeV/c^2)
  real(real64), parameter :: PROTON_MASS = 938.27208943_real64
  !> Neutron mass (MeV/c^2)
  real(real64), parameter :: NEUTRON_MASS = 939.56542194_real64
  !> Average nucleon mass (MeV/c^2)
  real(real64), parameter :: NUCLEON_MASS = (PROTON_MASS+NEUTRON_MASS)/2

  !> Reduced Planck constant (MeV·fm)
  real(real64), parameter :: HBAR = 197.3269804593025_real64

  !> Element names
  character(len=2), parameter :: ELEMENT_NAME(118) = ['H ','He','Li','Be','B ','C ', &
    & 'N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc', &
    & 'Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
    & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb', &
    & 'Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy', &
    & 'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl', &
    & 'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm', &
    & 'Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg', &
    & 'Cn','Nh','Fl','Mc','Lv','Ts','Og']


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Constants related to the dimension of the space and the quantum numbers
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !> Spatial dimension
  integer, parameter :: VDIM = 3 
  !> Index for directions
  integer, parameter :: VX = 1, VY = 2, VZ = 3
  !> Vector characters
  character(len=1), parameter :: VECTOR_CHAR(VDIM) = ['x','y','z']

  !> Spinor dimension (up or down)
  integer, parameter :: SDIM = 2 
  !> Index for spin up
  integer, parameter :: SU = 1, SD = 2
  !> Spin strings
  character(len=4), parameter :: SPIN_STR(SDIM) = ['up  ','down']
  !> Spin characters
  character(len=1), parameter :: SPIN_CHAR(SDIM) = ['u','d']
      
  !> Isospin dimension (proton or neutron)
  integer, parameter :: IDIM = 2
  !> Proton, neutron and matter dimension
  integer, parameter :: ADIM = IDIM + 1 
  !> Index for proton
  integer, parameter :: ISP = 1, ISN = 2, ISA = 3
  !> Isospin strings
  character(len=7), parameter :: ISOSPIN_STR(ADIM) = ['proton ','neutron','mass   ']
  !> Isospin characters
  character(len=1), parameter :: ISOSPIN_CHAR(ADIM) = ['p','n','a']

  !> Parity dimension (positive or negative)
  integer, parameter :: PDIM = 2
  !> Index for with and without parity operator
  integer, parameter :: OPI = 1, OPX = 2 ! OPI: no parity operator, OPX: parity operator
  !> Index for parity
  integer, parameter :: POS = 1, NEG = 2
  !> Parity strings
  character(len=3), parameter :: PARITY_STR(PDIM) = ['pos','neg']
  !> Parity characters
  character(len=1), parameter :: PARITY_CHAR(PDIM) =['+','-']

end module ModNuclConst
