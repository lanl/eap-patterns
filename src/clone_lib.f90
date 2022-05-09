! This is where MPI communications go, for now only serial
module clone_lib_module
  use define_kind
  use var_wrapper_class, only : var_wrapper
  public clone_get
  interface clone_get
     procedure clone_get_vw
     procedure clone_get_1
     procedure clone_get_2
  end interface clone_get
contains
  subroutine clone_get_vw(vw, intFlag)
    type(var_wrapper), intent(inout) :: vw(:)
    integer, intent(in) :: intFlag
  end subroutine clone_get_vw
  subroutine clone_get_1(invalue, value_work, nClone, nvars)
    real(REAL64), intent(in) :: invalue(:,:)
    real(REAL64), intent(out) :: value_work(:,:)
    integer, intent(in) :: nClone
    integer, intent(in) :: nvars
  end subroutine clone_get_1
  subroutine clone_get_2(myValue, nClone, nvars)
    real(REAL64), intent(inout) :: myValue(:,:)
    integer, intent(in) :: nClone
    integer, intent(in) :: nvars
  end subroutine clone_get_2
end module clone_lib_module

