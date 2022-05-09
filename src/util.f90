! ======================= copyright begin ======================== 
! Copyright (C) 2010-2013 Los Alamos National Security, LLC.            
! All rights Reserved.  See Copyright Notice File.                 
! Export Controlled Information                                    
! ======================== copyright end =========================
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

