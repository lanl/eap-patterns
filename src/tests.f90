module tests
  ! Tests of cell based loops

  public
  
contains
  function topcell_sum(m, values, n_iter) result(my_result)
    use iso_fortran_env, only: REAL64
    use mesh_types, only: mesh_t
    use clone_lib_module, only: clone_reduce, CLONE_SUM
    implicit none

    real(REAL64) :: my_result
    real(REAL64) :: partial_result
    
    type(mesh_t), intent(in) :: m
    real(REAL64), intent(in) :: values(:)
    integer, intent(in) :: n_iter

    integer :: i, iTop, iCell
    real(REAL64) :: local_sum

    my_result = 0.0_REAL64
    do i = 1, n_iter
       local_sum = 0.0_REAL64
       do iTop = 1, m%levels%numtop
          iCell = m%levels%ltop(iTop)
          local_sum = local_sum + values(iTop)
       end do
       call clone_reduce(partial_result, local_sum, CLONE_SUM)
       my_result = my_result + partial_result
    end do
       
  end function topcell_sum
end module tests
    
