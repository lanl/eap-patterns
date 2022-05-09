! ======================= copyright begin ========================       
! Copyright (C) 1990-2010 Los Alamos National Security, LLC.             
! All rights Reserved.  See Copyright Notice File.                       
!              and                                                       
! Copyright (C) 1990-2007 Science Applications International Corporation 
! Export Controlled Information                                          
! ======================== copyright end =========================       
! ------------------------------------------------------------------------------
module define_kind
  !*******************************************************************************
  !                                                                              *
  ! Defines generic parameters for data sizes, direction and side labels,        *
  ! useful numerical constants, and parallel processor identification.           *
  !                                                                              *
  !*******************************************************************************
  use iso_fortran_env, only : OUTPUT_UNIT, ERROR_UNIT, REAL64, REAL32, INT64, INT32, INT16
  !use, intrinsic :: ieee_arithmetic    !this has various query functionality for ieee_* conformance of fp variables 
  !*******************************************************************************
  ! We should move to the use of the standard definitions included in the module
  ! iso_fortran_env. However there are some name clashes that need to be resolved
  ! first, for instance in iso_fortran_env INT8 is an 8 bit rather than an 8 byte
  ! integer. In out usage we would want INT8. These types of global substitutions
  ! will need to be made before we can unreservedly use iso_fortran_env.
  !*******************************************************************************
  implicit none

  public

  integer, parameter :: INT4   = INT32
  integer, parameter :: INT8   = INT64

  integer, parameter :: REAL4  = REAL32
  integer, parameter :: REAL8  = REAL64

  integer, parameter :: CMPLX4  = SELECTED_REAL_KIND(6)
  integer, parameter :: CMPLX8  = SELECTED_REAL_KIND(12)

  integer, parameter :: CMPLX32 = SELECTED_REAL_KIND(6)
  integer, parameter :: CMPLX64 = SELECTED_REAL_KIND(12)

  integer, parameter :: WORDSIZE         = 8
  integer, parameter :: CHARS_PER_REAL64 = 8     !*almost* certainly always true.
  !
  ! Limits on filename length and full path length.
  ! Generally these are defined in limits.h (on linux, this may be
  ! found in /usr/include/linux/limits.h for example, but in similar 
  ! places on other OSs too) or more likely a machine dependent file
  ! included in limits.h.
  !
  integer, parameter :: LEN_PRBNM = 256 ! Changed from 32 to 256 to accommodate longer problem names
  ! Since strings of length LEN_PRBNM are
  ! printed to dumps this value should be a
  ! multiple of 16
  integer, parameter :: NAME_MAX      = 255+LEN_PRBNM
  integer, parameter :: PATH_MAX      = 4096+LEN_PRBNM
  integer, parameter :: HOST_NAME_MAX = 64         !see note in man gethostname(2)
  integer, parameter :: RZAID_MAX     = 32  ! longest-allowed reaction ZAID string; longest I found
  ! is 20-25 but round up to next power of two just in case

  real(REAL64), parameter :: HOUR =3600.0_REAL64

  integer, parameter :: ALL_DIR  = 0
  integer, parameter :: X_DIR    = 1
  integer, parameter :: Y_DIR    = 2
  integer, parameter :: Z_DIR    = 3

  integer, parameter :: LO_SIDE  = 1
  integer, parameter :: HI_SIDE  = 2

  integer, parameter :: LEN_LABS  = 32
  integer, parameter :: LEN_LABL  = 32
  integer, parameter :: LEN_UNITS = 32
  integer, parameter :: LEN_DANDT = 16
  integer, parameter :: LEN_DANDT_FULL = 20 ! with additional precision (ms)

  integer, parameter :: NOTHING_TO_DO = 0
  integer, parameter :: CRAY_TO_IEEE  = 1

  integer, parameter :: MAX_DEBUG = 50

  ! ---------------------
  ! HUGE and TINY flavors
  ! ---------------------
  ! EAP_HUGE_<type> and EAP_TINY_<type>
  ! EAP_HUGE = 1/EAP_TINY
  !
  ! These are defined to be reasonably large/small and be a reasonable
  ! distance away from HUGE(ONE), TINY(ONE).
  ! The goal is to be able to do some operations on these numbers and
  ! not get a FPE (floating point exception).
  ! We have (had?) code that did the following:
  !    val = HUGE(ONE)  (1.blahE+308 or something close to this)
  !    ...
  !    val_new = val * weight
  !    ...
  !    val = val_new / weight
  ! With finite precision math, sometimes val>HUGE(ONE) and we get FPE.
  ! So, best to use the below flavors because you never know who
  ! might take your value and try to do some operations on it.
  ! 
  ! If you think you really to have a REAL32 flavor...think instead
  ! about just converting your data type to REAL64.  Sure, there might
  ! be some performance/memory reason for it...but we are a double
  ! precision code.
  real(REAL64), parameter :: EAP_HUGE_REAL64 = 1.0d+250
  real(REAL64), parameter :: EAP_TINY_REAL64 = 1.0d-250

  ! Numbers REAL64
  real(REAL64), parameter :: ZERO      =    0.0_REAL64
  real(REAL64), parameter :: EIGHTH    =    0.125_REAL64
  real(REAL64), parameter :: QUARTER   =    0.25_REAL64
  real(REAL64), parameter :: HALF      =    0.5_REAL64
  real(REAL64), parameter :: ONE       =    1.0_REAL64
  real(REAL64), parameter :: TWO       =    2.0_REAL64
  real(REAL64), parameter :: THREE     =    3.0_REAL64
  real(REAL64), parameter :: FOUR      =    4.0_REAL64
  real(REAL64), parameter :: FIVE      =    5.0_REAL64
  real(REAL64), parameter :: SIX       =    6.0_REAL64
  real(REAL64), parameter :: SEVEN     =    7.0_REAL64
  real(REAL64), parameter :: EIGHT     =    8.0_REAL64
  real(REAL64), parameter :: NINE      =    9.0_REAL64
  real(REAL64), parameter :: TEN       =   10.0_REAL64
  real(REAL64), parameter :: ELEVEN    =   11.0_REAL64
  real(REAL64), parameter :: TWELVE    =   12.0_REAL64
  real(REAL64), parameter :: FIFTEEN   =   15.0_REAL64
  real(REAL64), parameter :: TWENTY    =   20.0_REAL64
  real(REAL64), parameter :: HUNDRED   =  100.0_REAL64
  real(REAL64), parameter :: THOUSAND  = 1000.0_REAL64
  real(REAL64), parameter :: MILLION   = 1.0e6_REAL64
  real(REAL64), parameter :: NEGONE    = -1.0_REAL64

  real(REAL64), parameter :: EGAMMA    = 0.57721566490153286061_REAL64 ! Euler's gamma
  real(REAL64), parameter :: PI        = 3.14159265358979323846_REAL64
  real(REAL64), parameter :: TWOPI     = TWO*PI
  real(REAL64), parameter :: FOURPI    = FOUR*PI
  real(REAL64), parameter :: FOURPIO3  = FOURPI/THREE
  real(REAL64), parameter :: THIRD     = ONE/THREE
  real(REAL64), parameter :: TWOTHIRD  = TWO/THREE
  real(REAL64), parameter :: FOURTHIRD = FOUR/THREE
  real(REAL64), parameter :: FIFTH     = ONE/FIVE
  real(REAL64), parameter :: SIXTH     = ONE/SIX
  real(REAL64), parameter :: TENTH     = ONE/TEN
  real(REAL64), parameter :: THREEHALF  = THREE/TWO
  real(REAL64), parameter :: FIVEFOURTH = FIVE/FOUR

  ! Numbers INT64
  integer(INT64),  parameter :: ZEROI8     = 0_INT64
  integer(INT64),  parameter :: ONEI8      = 1_INT64
  integer(INT64),  parameter :: TWOI8      = 2_INT64
  integer(INT64),  parameter :: THREEI8    = 3_INT64
  integer(INT64),  parameter :: FOURI8     = 4_INT64
  integer(INT64),  parameter :: FIVEI8     = 5_INT64
  integer(INT64),  parameter :: SIXI8      = 6_INT64
  integer(INT64),  parameter :: SEVENI8    = 7_INT64
  integer(INT64),  parameter :: EIGHTI8    = 8_INT64
  integer(INT64),  parameter :: NEGONEI8   = -1_INT64
  integer(INT64),  parameter :: THOUSANDI8 = 1000_INT64

  integer    , save :: io_stat ,    &    ! status returned by ciolib functions != 0 on fail
       mype    ,    &    ! identifer for local processor number
       numpe   ,    &    ! number of processors being used
       maxpe   ,    &    ! number of last processor = numpe-1
       iope    ,    &    ! identifier of processor performing IO.
       mpierror,    &    ! used to test for mpi errors
       allostat          ! used to test allocate status

  integer     , save :: define_kind_debug = 0
  integer     , save :: define_kind_i (MAX_DEBUG) = 0
  real(REAL64), save :: define_kind_r8(MAX_DEBUG) = ZERO


  ! version information
  ! version info
  character(LEN=:), allocatable :: &
       v_special_version, &
       v_compile_date, &
       v_tf_prebuilt_name, &
       scm_id_0, &
       scm_id_1, &
       scm_id_2, &
       scm_id_3, &
       scm_id_4, &
       l_eap_version

  ! ==============================================================================
contains
  ! ------------------------------------------------------------------------------
  subroutine test_allostat(comment, line, file)
    !*******************************************************************************
    !                                                                              *
    ! Test allostat and abort in not zero                                          *
    !                                                                              *
    !*******************************************************************************
    character(*), intent(in), optional :: comment
    integer, intent(in), optional :: line
    character(*), intent(in), optional :: file


    if(allostat.ne.0)then
       stop 'unable to allocate'
    endif

  end subroutine test_allostat

  ! ============================================================================

  ! ------------------------------------------------------------------------------
end module define_kind
! ==============================================================================
