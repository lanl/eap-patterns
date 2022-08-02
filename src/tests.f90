module tests
  ! Tests of cell based loops
  use iso_fortran_env, only: INT64, REAL64
  use clone_lib_module, only: clone_barrier, clone_reduce, CLONE_SUM
  use pio_interface, only: pio_now
  public
  integer, private :: myid
contains
  subroutine test_driver(m, n_iter)
    use mesh_types, only: mesh_t
    use clone_lib_module, only: clone_myid
    implicit none
    type(mesh_t), intent(in) :: m
    integer, intent(in) :: n_iter

    myid = clone_myid()
    
    if (myid == 0 ) write(*,'(/,"-------BEGIN TESTS----------",/)')
    call topcell_sum(m, n_iter)
    call faces_sum(m, n_iter)
    if (myid == 0 ) write(*,'(/,"--------END TESTS-----------",/)')
  end subroutine test_driver

  subroutine printit(the_name, the_status, the_dt)
    use clone_lib_module, only: clone_myid
    implicit none
    character(len=*), intent(in) :: the_name
    logical, intent(in) :: the_status
    real(REAL64) :: the_dt
    if (clone_myid() == 0) then
       if (the_status) then
          write(*,*) '    PASS: ', the_dt, trim(the_name)
       else
          write(*,*) '  **FAIL: ', the_dt, trim(the_name)
       end if
    end if
  end subroutine printit

  subroutine topcell_sum(m, n_iter) 
    use iso_fortran_env, only: REAL64
    use mesh_types, only: mesh_t
    implicit none

    type(mesh_t), intent(in) :: m
    
    real(REAL64), allocatable :: values(:)
    integer, intent(in) :: n_iter

    integer :: i, iTop, iCell
    real(REAL64) :: my_dt, local_sum, my_sum, partial_result, expected_result

    ! Initialize arrays
    allocate(values(m%cells%numcell))
    values = 1
    
    ! Compute the expected result
    partial_result = real(m%levels%numtop, kind=REAL64)
    call clone_reduce(expected_result, partial_result, CLONE_SUM)
    expected_result = real(n_iter,kind=REAL64) * expected_result

    ! Run the loop
    my_sum = 0.0_REAL64
    call clone_barrier()
    my_dt = pio_now() 
    call clone_barrier()
    do i = 1, n_iter
       local_sum = 0.0_REAL64
       partial_result = 0.0_REAL64
       do iTop = 1, m%levels%numtop
          iCell = m%levels%ltop(iTop)
          local_sum = local_sum + values(iTop)
       end do
       call clone_reduce(partial_result, local_sum, CLONE_SUM)
       my_sum = my_sum + partial_result
    end do
    call clone_barrier()
    my_dt = pio_now() - my_dt

    ! Report the results
    if (myid == 0) then
       if (my_sum /= expected_result) then
          write(*,*) "    topcell_sum fail: ",my_sum, " != ", expected_result
       end if
    end if
    call printit("topcell_sum", (my_sum == expected_result), my_dt)
  end subroutine topcell_sum
  
  subroutine faces_sum(m, n_iter)
    use define_kind, only: HI_SIDE, LO_SIDE
    use iso_fortran_env, only: REAL64, INT64
    use mesh_types, only: mesh_t
    implicit none

    type(mesh_t), intent(in) :: m
    
    real(REAL64), allocatable :: values(:)
    integer, intent(in) :: n_iter

    integer :: i, iTop, iCell, iFace, iDim, iLoop, iType, n, nlo, nhi
    real(REAL64) :: my_dt
    integer :: faces_by_types(5), faces_on_pe_boundary(5)
    integer(INT64) :: global_faces_by_types(5), global_faces_on_pe_boundary(5)

    ! Initialize arrays
    allocate(values(m%cells%numcell))
    values = 1

    faces_by_types = 0
    faces_on_pe_boundary = 0
    global_faces_by_types = 0
    call clone_barrier()
    my_dt = pio_now() 

    do iDim = 1, m%sim%numdim
       ! for each dimension
       do iLoop = 1, m%faces%face_num(iDim)
          ! for each face type in that dimension
          iType = m%faces%face_id(iLoop, iDim)
          if (iType <= 2 ) then
             ! Types 1 & 2 are always local
             n = 1 + m%faces%face_hi(iLoop, iDim) - m%faces%face_lo(iLoop, iDim)
             faces_by_types(iType) = faces_by_types(iType) + n
          else if ( iType > 2 ) then
             ! Check for PE boundary faces
             nlo = m%faces%face_lo(iLoop, iDim)
             nhi = m%faces%face_hi(iLoop, iDim)
             do n = nlo, nhi
                ! if (m%faces%face_local(n, LO_SIDE, idim) > m%cells%numcell) then
                !    faces_on_pe_boundary(iType) = faces_on_pe_boundary(iType) + 1
                ! end if
                ! All low side faces belong to us
                faces_by_types(iType) = faces_by_types(iType) + 1
                
                if (m%faces%face_local(n, HI_SIDE, idim) > m%cells%numcell) then
                   faces_on_pe_boundary(iType) = faces_on_pe_boundary(iType) + 1
                else
                   ! Only non-boundary high side faces belong to us
                   faces_by_types(iType) = faces_by_types(iType) + 1
                end if
             end do
          end if
       end do
    end do

    do iType = 1, 5
       ! Compute the number of faces
       call clone_reduce(global_faces_by_types(iType), faces_by_types(iType), CLONE_SUM)
       call clone_reduce(global_faces_on_pe_boundary(iType), faces_on_pe_boundary(iType), CLONE_SUM)
    end do
    call clone_barrier()
    my_dt = pio_now() - my_dt

    ! Report the results
    if (myid == 0) then
       write(*,*) 'type   real_faces pe_boundary_faces'
       do iType = 1, 5
          write(*,'(I4, " ", 2(I12, " "))') iType, global_faces_by_types(iType)- global_faces_on_pe_boundary(iType), &
               global_faces_on_pe_boundary(iType)
       end do
    end if
    call printit("faces_sum", .true., my_dt)
  end subroutine faces_sum
end module tests
    
