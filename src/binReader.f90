!!  =======================================================================
!!  © (or copyright) 2022. Triad National Security, LLC. All rights
!!  reserved.  This program was produced under U.S. Government contract
!!  89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
!!  operated by Triad National Security, LLC for the U.S.  Department of
!!  Energy/National Nuclear Security Administration. All rights in the
!!  program are reserved by Triad National Security, LLC, and the
!!  U.S. Department of Energy/National Nuclear Security
!!  Administration. The Government is granted for itself and others acting
!!  on its behalf a nonexclusive, paid-up, irrevocable worldwide license
!!  in this material to reproduce, prepare derivative works, distribute
!!  copies to the public, perform publicly and display publicly, and to
!!  permit others to do so.
!!
!!  See LICENSE file for details
!!  =======================================================================

#define C0(x) (trim(x)//char(0))

module binreader
  use iso_fortran_env, only: REAL64, INT64, INT32, INT8
  use iso_c_binding
  implicit none

  private

  public :: binfile
  public :: binfile_verify_signature
  
  type :: variable_t
     character(len=128) :: name
     integer(INT64) :: size
     integer(INT64) :: offset
  end type variable_t

  type :: binfile
     integer(INT64) :: ndim
     integer(INT64) :: nMat
     integer(INT64) :: nCells
     integer(INT64) :: nVars
     type(variable_t), pointer :: vars(:)
     integer(c_int64_t) :: fp
   contains
     procedure :: init => binfile_init
     procedure :: binfile_read_i8_2
     procedure :: binfile_read_i64_2
     procedure :: binfile_read_f64_2
     procedure :: binfile_read_i32
     procedure :: binfile_read_i64
     procedure :: binfile_read_f64
     generic :: read => binfile_read_i8_2,  &
                        binfile_read_i64_2, &
                        binfile_read_f64_2, &
                        binfile_read_i32,   &
                        binfile_read_i64,   &
                        binfile_read_f64
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
    
    call readIt(c_loc(iTmp), this%fp, offset, 8_INT64)
    this%ndim = iTmp
    offset = offset + 8
    
    call readIt(c_loc(iTmp), this%fp, offset, 8_INT64)
    this%nCells = iTmp
    offset = offset + 8
    
    call readIt(c_loc(iTmp), this%fp, offset, 8_INT64)
    this%nMat = iTmp
    offset = offset + 8
    
    call readIt(c_loc(vLen), this%fp, offset, 8_INT64)
    offset = offset + 8
    
    call readIt(c_loc(iTmp), this%fp, offset, 8_INT64)
    this%nVars = iTmp
    offset = offset + 8

    allocate(this%vars(this%nVars))

    BLOCK
      character(len=vLen), target :: tmpname
      integer(INT64), target :: tmpOffset, tmpSize
      integer :: ivar
      
      do ivar = 1, this%nVars
         call readIt(c_loc(tmpName), this%fp, offset, vLen)
         offset = offset + vLen
         call readIt(c_loc(tmpSize), this%fp, offset, 8_INT64)
         offset = offset + 8
         call readIt(c_loc(tmpOffset), this%fp, offset, 8_INT64)
         offset = offset + 8
         this%vars(ivar)%name = tmpname
         this%vars(ivar)%size = tmpSize
         this%vars(ivar)%offset = tmpOffset
      end do
    END BLOCK

  END subroutine binfile_init

  subroutine binfile_read_i8_2(this, outPtr, var, iStart, nCount, n2)
    ! different from other reads in that it reads i8 variable
    implicit none
    class(binfile) :: this
    integer(INT8), pointer, dimension(:,:), intent(inout) :: outPtr
    character(len=*), intent(in) :: var
    integer(INT64), intent(in) :: iStart, nCount, n2
    integer(INT64), target :: myStart, myBytes, myN
    integer :: i

    do i =1, this%nvars
       if (var == this%vars(i)%name) then
          myN = nCount
          myStart = this%vars(i)%offset + (iStart - 1) * n2 * this%vars(i)%size
          myBytes = myN * n2 * this%vars(i)%size
          if (.not. associated(outPtr)) then
             allocate(outPtr(n2,myN))
          end if
          call readIt(c_loc(outPtr), this%fp, myStart, myBytes)
          exit
       end if
    end do
  end subroutine binfile_read_i8_2

  subroutine binfile_read_i64_2(this, outPtr, var, iStart, nCount, n2)
    ! different from other reads in that it reads i8 variable
    implicit none
    class(binfile) :: this
    integer(INT64), pointer, dimension(:,:), intent(inout) :: outPtr
    character(len=*), intent(in) :: var
    integer(INT64), intent(in) :: iStart, nCount, n2
    integer(INT64), target :: myStart, myBytes, myN
    integer :: i
    do i =1, this%nvars
       if (var == this%vars(i)%name) then
          myN = nCount
          myStart = this%vars(i)%offset + (iStart - 1) * n2 * this%vars(i)%size
          myBytes = myN * n2 * this%vars(i)%size
          if (.not. associated(outPtr)) then
             allocate(outPtr(n2,myN))
          end if
          call readIt(c_loc(outPtr), this%fp, myStart, myBytes)
          exit
       end if
    end do
  end subroutine binfile_read_i64_2

  subroutine binfile_read_f64_2(this, outPtr, var, iStart, nCount, n2)
    ! different from other reads in that it reads i8 variable
    implicit none
    class(binfile) :: this
    real(REAL64), pointer, dimension(:,:), intent(inout) :: outPtr
    character(len=*), intent(in) :: var
    integer(INT64), intent(in) :: iStart, nCount, n2
    integer(INT64), target :: myStart, myBytes, myN
    integer :: i

    do i =1, this%nvars
       if (var == this%vars(i)%name) then
          myN = nCount
          myStart = this%vars(i)%offset + (iStart - 1) * n2 * this%vars(i)%size
          myBytes = myN * n2 * this%vars(i)%size
          if (.not. associated(outPtr)) then
             allocate(outPtr(n2,myN))
          end if
          call readIt(c_loc(outPtr), this%fp, myStart, myBytes)
          exit
       end if
    end do
  end subroutine binfile_read_f64_2

  subroutine binfile_read_i32(this, outPtr, var, iStart, nCount)
    ! different from other reads in that it reads i8 variable
    implicit none
    class(binfile) :: this
    integer(INT32), pointer, dimension(:), intent(inout) :: outPtr
    character(len=*), intent(in) :: var
    integer(INT64), intent(in) :: iStart, nCount
    integer(INT64), target :: myStart, myBytes, myN
    integer :: i

    do i =1, this%nvars
       if (var == this%vars(i)%name) then
          myN = nCount
          myStart = this%vars(i)%offset + (iStart - 1) * this%vars(i)%size
          myBytes = myN * this%vars(i)%size
          if (.not. associated(outPtr)) then
             allocate(outPtr(myN))
          end if
          call readIt(c_loc(outPtr), this%fp, myStart, myBytes)
          exit
       end if
    end do
  end subroutine binfile_read_i32

  subroutine binfile_read_i64(this, outPtr, var, iStart, nCount)
    ! different from other reads in that it reads i8 variable
    implicit none
    class(binfile) :: this
    integer(INT64), pointer, dimension(:), intent(inout) :: outPtr
    character(len=*), intent(in) :: var
    integer(INT64), intent(in) :: iStart, nCount
    integer(INT64), target :: myStart, myBytes, myN
    integer :: i

    do i =1, this%nvars
       if (var == this%vars(i)%name) then
          myN = nCount
          myStart = this%vars(i)%offset + (iStart - 1) * this%vars(i)%size
          myBytes = myN * this%vars(i)%size
          if (.not. associated(outPtr)) then
             allocate(outPtr(myN))
          end if
          call readIt(c_loc(outPtr), this%fp, myStart, myBytes)
          return
       end if
    end do
    stop 'variable not found'
  end subroutine binfile_read_i64

  subroutine binfile_read_f64(this, outPtr, var, iStart, nCount)
    ! different from other reads in that it reads i8 variable
    implicit none
    class(binfile) :: this
    real(REAL64), pointer, dimension(:), intent(inout) :: outPtr
    character(len=*), intent(in) :: var
    integer(INT64), intent(in) :: iStart, nCount
    integer(INT64), target :: myStart, myBytes, myN
    integer :: i

    do i =1, this%nvars
       if (var == this%vars(i)%name) then
          myN = nCount
          myStart = this%vars(i)%offset + (iStart - 1) * this%vars(i)%size
          myBytes = myN * this%vars(i)%size
          if (.not. associated(outPtr)) then
             allocate(outPtr(myN))
          end if
          call readIt(c_loc(outPtr), this%fp, myStart, myBytes)
          exit
       end if
    end do
  end subroutine binfile_read_f64
 
end module binreader

