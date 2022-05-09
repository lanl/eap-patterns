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

module util
  implicit none
  public
contains
  subroutine global_error(s)
    !*******************************************************************************
    !                                                                              *
    ! report location of error and abort when ALL PEs have the same error          *
    !                                                                              *
    !*******************************************************************************
    character*(*), intent(in) :: s
    write(*,610)s
    stop 'stop'

610 format('GLOBAL_ERROR called: ',a)

  end subroutine global_error

end module util

