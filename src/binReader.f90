#define C0(x) (trim(x)//char(0))

module binreader
  use iso_fortran_env, only: REAL64, INT64, INT32
  use iso_c_binding
  implicit none

  private

  public :: binfile
  public :: binfile_verify_signature
  
  type :: variable_t
     character(len=128) :: name
     integer(INT64) :: offset
  end type variable_t

  type :: binfile
     integer(INT64) :: ndim
     integer(INT64) :: nCells
     integer(INT64) :: nVars
     type(variable_t), pointer :: vars(:)
     integer(c_int64_t) :: fp
   contains
     procedure :: init => binfile_init
     procedure :: read_i64 => binfile_read_i64
     procedure :: read_i64_2d => binfile_read_i64_2
     procedure :: read_f64 => binfile_read_f64
     procedure :: release => binfile_release
  end type binfile

  ! Interface to C subroutines
  interface
     function openIt(fname) BIND(C, name="openIt")
       use iso_c_binding
       implicit none
       integer(c_int64_t) :: openIt
       character(c_char), intent(in) :: fname(1)
     end function openIt

     subroutine seekIt(fp, pos) BIND(C, name="seekIt")
       use iso_c_binding
       implicit none
       integer(c_int64_t), VALUE :: fp
       integer(c_int64_t) :: pos
     end subroutine seekIt

     subroutine closeIt(fp) BIND(C, name="closeIt")
       use iso_c_binding
       implicit none
       integer(c_int64_t), VALUE :: fp
     end subroutine closeIt

     subroutine readIt(ptr, fp, offset, isize) BIND(C, name="readIt")
       use iso_c_binding
       implicit none
       type(c_ptr), VALUE :: ptr
       integer(c_int64_t), VALUE  :: fp
       integer(c_int64_t), VALUE, intent(in) :: offset
       integer(c_int64_t), VALUE, intent(in) :: isize
     end subroutine readIt
  end interface
  
contains

  subroutine binfile_release(this)
    implicit none
    class(binfile) :: this
    call closeIt(this%fp)
  end subroutine binfile_release
  
  logical function binfile_verify_signature(fname)
    use iso_c_binding
    implicit none
    character(len=*) :: fname
    character(len=16), target :: sig
    integer(INT64) :: fp
    fp = openIt(C0(fname))
    call readIt(c_loc(sig), fp, 0_INT64, 16_INT64)
    call closeIt(fp)
    binfile_verify_signature = ( sig == 'eap-patterns-bin')
  end function binfile_verify_signature
  
  subroutine binfile_init(this, fname)
    implicit none
    class(binfile) :: this
    character(len=*), intent(in) :: fname
    character(16), target :: sig
    integer(INT64), target :: iEndian, vLen, iTmp, offset


    this%fp = openIt(C0(fname))

    offset = 0
    call readIt(c_loc(sig), this%fp, offset, 16_INT64)
    offset = offset + 16
    call readIt(c_loc(iEndian), this%fp, offset, 8_INT64)
    offset = offset + 8
    call readIt(c_loc(iTmp), this%fp, 24_INT64, 8_INT64)
    offset = offset + 8
    this%ndim = iTmp
    call readIt(c_loc(iTmp), this%fp, 32_INT64, 8_INT64)
    offset = offset + 8
    this%nCells = iTmp
    call readIt(c_loc(vLen), this%fp, 40_INT64, 8_INT64)
    offset = offset + 8
    call readIt(c_loc(iTmp), this%fp, 48_INT64, 8_INT64)
    offset = offset + 8
    this%nVars = iTmp

    allocate(this%vars(this%nVars))

    BLOCK
      character(len=vLen), target :: tmpname
      integer(INT64), target :: tmpoffset
      integer :: ivar
      
      do ivar = 1, this%nVars
         call readIt(c_loc(tmpName), this%fp, offset, vLen)
         offset = offset + vLen
         call readIt(c_loc(tmpoffset), this%fp, offset, 8_INT64)
         offset = offset + 8
         this%vars(ivar)%name = tmpname
         this%vars(ivar)%offset = tmpoffset
      end do
    END BLOCK

  END subroutine binfile_init

  function binfile_read_i64(this, var, iStart, nCount) result(read_i64)
    implicit none
    real(c_double), pointer, dimension(:) :: tmpPtr
    integer(INT64), pointer, dimension(:) :: read_i64
    class(binfile) :: this
    character(len=*), intent(in) :: var
    integer(INT64), intent(in) :: iStart, nCount
    integer(INT64), target :: myStart, myN, i

    nullify(read_i64)
    nullify(tmpPtr)

    do i =1, this%nvars
       if (var == this%vars(i)%name) then
          myStart = iStart
          myN = nCount
          allocate(tmpPtr(myN))
          allocate(read_i64(myN))
          call readIt(c_loc(tmpPtr), this%fp, this%vars(i)%offset + (iStart - 1) * 8, 8_INT64 * nCount)
          exit
       end if
    end do
    if (associated(read_i64)) then
       do i = 1, nCount
          read_i64(i) = int(tmpPtr(i), kind=INT64)
       end do
    end if

    if (associated(tmpPtr)) then
       deallocate(tmpPtr)
    end if
       
    return
  end function binfile_read_i64
  
  function binfile_read_i64_2(this, var, n2, iStart, nCount) result(read_i64_2)
    implicit none
    real(c_double), pointer, dimension(:) :: tmp
    integer(INT64), pointer, dimension(:,:) :: read_i64_2
    integer(INT64), intent(in) :: iStart, nCount
    class(binfile) :: this
    character(len=*), intent(in) :: var
    integer, intent(in) :: n2
    character(len=32) :: suffix
    character(len=128) :: tmpChar
    integer :: iVar, j, iCell
    integer(INT64), target :: myN, myStart
    integer(INT64) :: offset


    myStart = iStart
    myN = nCount
    allocate(tmp(myN))
    allocate(read_i64_2(myN, n2))
    do j = 1, n2
       write(suffix,*) j
       tmpChar  = trim(var) // '_' // trim(adjustl(suffix))
       VARLOOP: do iVar =1, this%nvars
          if (trim(tmpChar) == this%vars(iVar)%name) then
             offset = this%vars(iVar)%offset + (iStart - 1) * 8
             call readIt(c_loc(tmp), this%fp, offset, 8_INT64 * nCount)
             do iCell = 1, nCount
                read_i64_2(iCell, j) = int(tmp(iCell), kind=INT64)
             end do
             exit VARLOOP
          end if
       end do VARLOOP
    end do

    deallocate(tmp)
       
    return
  end function binfile_read_i64_2
 
  function binfile_read_f64(this, var, iStart, nCount) result(read_f64)
    implicit none
    real(REAL64), pointer, dimension(:) :: read_f64
    class(binfile) :: this
    character(len=*), intent(in) :: var
    integer(INT64), intent(in) :: iStart, nCount
    integer(INT64), target :: myStart, myN, i

    nullify(read_f64)

    do i =1, this%nvars
       if (var == this%vars(i)%name) then
          myStart = iStart
          myN = nCount
          allocate(read_f64(myN))
          call readIt(c_loc(read_f64), this%fp, this%vars(i)%offset + (iStart - 1) * 8, 8_INT64 * nCount)
          exit
       end if
    end do
    return
  end function binfile_read_f64
 
end module binreader

