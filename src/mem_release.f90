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

! a simple interface to release memory
module mem_release
  use iso_fortran_env, only : REAL64, INT64
  use iso_c_binding, only: c_bool
  implicit none
  public
  interface release
     module procedure release_d_1
     module procedure release_d_2

     module procedure release_i
     module procedure release_i_1
     module procedure release_i_2
     module procedure release_i_3

     module procedure release_i64
     module procedure release_i64_1
     module procedure release_i64_2
     module procedure release_i64_3

     module procedure release_l
     module procedure release_l_1
  end interface release
contains

  !---------Doubles------------
  subroutine release_d(ptr)
    REAL(REAL64), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_d

  subroutine release_d_1(ptr)
    REAL(REAL64), dimension(:), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_d_1

  subroutine release_d_2(ptr)
    REAL(REAL64), dimension(:,:), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_d_2

  !---------Integer 32------------ 
  subroutine release_i(ptr)
    integer, pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_i

  subroutine release_i_1(ptr)
    integer, dimension(:), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_i_1

  subroutine release_i_2(ptr)
    integer, dimension(:,:), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_i_2

  subroutine release_i_3(ptr)
    integer, dimension(:,:,:), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_i_3

  !---------Integer 64------------
  subroutine release_i64(ptr)
    integer(INT64), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_i64

  subroutine release_i64_1(ptr)
    integer(INT64), dimension(:), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_i64_1

  subroutine release_i64_2(ptr)
    integer(INT64), dimension(:,:), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_i64_2

  subroutine release_i64_3(ptr)
    integer(INT64), dimension(:,:,:), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_i64_3

  !---------Logicals------------
  subroutine release_l(ptr)
    logical(c_bool), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_l

  subroutine release_l_1(ptr)
    logical(c_bool), dimension(:), pointer :: ptr
    if (associated(ptr)) deallocate(ptr)
    nullify(ptr)
  end subroutine release_l_1

end module mem_release

