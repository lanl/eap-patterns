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

module tests
  ! Tests of cell based loops
  use iso_fortran_env, only: INT64, REAL64
  use clone_lib_module, only: clone_barrier, clone_reduce, CLONE_SUM, CLONE_MAX, clone_myid
  public
  integer, private :: myid
contains

  real(REAL64) function now()
    integer(INT64) :: c_now
    real(REAL64) :: c_rate
    call system_clock(c_now, c_rate)
    now = real(c_now, REAL64) / c_rate
  end function now
    
  subroutine test_driver(fm, n_iter)
    use fakemesh, only: fakemesh_t
    use clone_lib_module, only: clone_myid
    implicit none
    type(fakemesh_t) :: fm
    integer, intent(in) :: n_iter

    myid = clone_myid()
    
    if (myid == 0 ) write(*,'(/,"-------BEGIN TESTS----------",/)')
    ! call faces_sum(fm%m, n_iter)
    ! call faces_scatter(fm%m, n_iter)
    ! call topcell_sum(fm%m, n_iter)

    call faces_check(fm%m, n_iter)

    ! The true test: derivatives
    call deriv_test(fm, n_iter)
    
    if (myid == 0 ) write(*,'(/,"--------END TESTS-----------",/)')
  end subroutine test_driver

  subroutine printit(the_name, the_status, the_dt)
    use clone_lib_module, only: clone_myid
    implicit none
    character(len=*), intent(in) :: the_name
    logical, intent(in) :: the_status
    real(REAL64) :: the_dt, max_dt

    call clone_reduce(max_dt, the_dt, CLONE_MAX)
    if (clone_myid() == 0) then
       if (the_status) then
          write(*,*) '    PASS: ', max_dt, trim(the_name)
       else
          write(*,*) '  **FAIL: ', max_dt, trim(the_name)
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
    expected_result = 1.0d0
    partial_result = real(m%levels%numtop, kind=REAL64)
    call clone_reduce(expected_result, partial_result, CLONE_SUM, .false.)
    expected_result = real(n_iter,kind=REAL64) * expected_result

    ! Run the loop
    my_sum = 0.0_REAL64
    call clone_barrier()
    my_dt = now() 
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
    my_dt = now() - my_dt

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
    real(REAL64) :: my_dt, delta, factor
    real(REAL64) :: faces_by_types(5), faces_on_pe_boundary(5)
    real(REAL64) :: global_faces_by_types(5), global_faces_on_pe_boundary(5)

    ! Initialize arrays
    allocate(values(m%cells%numcell))
    values = 1.0d0

    faces_by_types = 0
    faces_on_pe_boundary = 0
    global_faces_by_types = 0
    call clone_barrier()
    my_dt = now() 

    if (m%sim%numdim == 1) then
       factor = 1.0d0
    else if (m%sim%numdim == 2) then
       factor = 1.0d0
    else if (m%sim%numdim == 3) then
       factor = 1.25d0
    end if
    
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
             if ( iType == 3 ) then
                delta = 2.0d0
             else
                delta = factor
             end if
             do n = nlo, nhi
                ! only low side boundary faces belong to us
                if (m%faces%face_local(n, HI_SIDE, idim) > m%cells%numcell) then
                   faces_on_pe_boundary(iType) = faces_on_pe_boundary(iType) + delta
                   faces_by_types(iType) = faces_by_types(iType) + delta
                else
                   faces_by_types(iType) = faces_by_types(iType) + delta
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
    my_dt = now() - my_dt

    ! Report the results
    if (myid == 0) then
       do iType = 1, 5
          !write(*,'(I4, " ", 2(I12, " "))') iType, global_faces_by_types(iType)- global_faces_on_pe_boundary(iType), &
          ! global_faces_on_pe_boundary(iType)
       end do
       write(*,'("            Faces: ", 5(F12.1, " "))') &
         (global_faces_by_types(iType)- global_faces_on_pe_boundary(iType), iType=1,5)
       write(*,'("       Bdry Faces: ", 5(F12.1, " "))') &
         (global_faces_on_pe_boundary(iType), iType=1,5)
    end if
    call printit("faces_sum", .true., my_dt)
  end subroutine faces_sum

  subroutine faces_check(m, n_iter)
    use define_kind, only: HI_SIDE, LO_SIDE
    use iso_fortran_env, only: REAL64, INT64
    use mesh_types, only: mesh_t
#ifdef ENABLE_VTUNE  
    use ittnotify
#endif
    implicit none

    type(mesh_t), intent(in) :: m
    
    real(REAL64), allocatable :: values(:)
    integer, intent(in) :: n_iter

    integer :: i, iTop, iCell, iFace, iDim, iLoop, iType, n, nlo, nhi, ilo, ihi, iIter
    real(REAL64) :: my_dt, my_sum, expected_result, my_T_sum, all_sum, all_expected
    integer, allocatable :: mothers(:)
    

    ! Initialize arrays
    allocate(values(m%cells%numcell_clone))
    values = 0.0d0

#ifdef ENABLE_VTUNE  
    call itt_resume()
#endif
    call clone_barrier()
    my_dt = now()
    if (m%sim%numdim == 2) then
       my_T_sum = 0.5
    else if (m%sim%numdim == 3) then
       my_T_sum = 0.25
    else
       my_T_sum = 1.0
    end if
    do iIter = 1, n_iter
       do iDim = 1, m%sim%numdim
          ! for each dimension
          do iLoop = 1, m%faces%face_num(iDim)
             ! for each face type in that dimension
             iType = m%faces%face_id(iLoop, iDim)
             nlo = m%faces%face_lo(iLoop, iDim)
             nhi = m%faces%face_hi(iLoop, iDim)
             do n = nlo, nhi
                ilo = m%faces%face_local(n, LO_SIDE, idim)
                ihi = m%faces%face_local(n, HI_SIDE, idim)
                if (iType <= 2 ) then
                   ! Types 1, & 2 have same cell on both sides
                   if (ilo >= 1 .and. ilo <= m%cells%numcell) then
                      values(ilo) = values(ilo) + 1.0D0
                   else
                      write(*,*) clone_myid(), ': ilo killer:', n, iType, iDim, ihi, ilo
                   end if
                else if ( iType == 3 ) then
                   if (ilo <= m%cells%numcell) then
                      values(ilo) = values(ilo) + 1.0D0
                   end if
                   if (ihi <= m%cells%numcell) then
                      values(ihi) = values(ihi) + 1.0D0
                   end if
                else if ( iType == 4 ) then
                   if (ilo <= m%cells%numcell) then
                      values(ilo) = values(ilo) + 1.0D0
                   end if
                   if (ihi <= m%cells%numcell) then
                      values(ihi) = values(ihi) + my_T_sum
                   end if
                else if ( iType == 5 ) then
                   if (ilo <= m%cells%numcell) then
                      values(ilo) = values(ilo) + my_T_sum
                   end if
                   if (ihi <= m%cells%numcell) then
                      values(ihi) = values(ihi) + 1.0D0
                   end if
                end if
             end do
          end do
       end do
    end do
    call clone_barrier()
#ifdef ENABLE_VTUNE  
    call itt_pause()
#endif
    my_dt = now() - my_dt

    ! Check results for leaf cells
    allocate(mothers(m%cells%numcell_clone))
    mothers = 1.0
    my_sum = 0.0D0
    expected_result = real(n_iter * 2 * m%sim%numdim, REAL64) * m%levels%numtop
    do i = 1, m%levels%numtop
       iCell = m%levels%ltop(i)
       mothers(iCell) = 0.0
       my_sum = my_sum + values(iCell)
       ! reset values(i) to known bad value
    end do

    call clone_reduce(all_sum, my_sum, CLONE_SUM)
    call clone_reduce(all_expected, expected_result, CLONE_SUM)
    call clone_barrier()
    if (clone_myid() == 0 .and. all_expected /= all_sum) then
       write(*,*) clone_myid(), "daughters:", all_sum, all_expected, all_sum - all_expected
    end if
    call printit("face_check daughters", (all_sum ==  all_expected), my_dt)

    ! check results for mother cells
    my_sum = 0.0
    do iCell = 1, m%cells%numcell
       my_sum = my_sum + mothers(iCell) * values(iCell)
    end do
    expected_result = 0.0D0
    call clone_reduce(all_sum, my_sum, CLONE_SUM)
    call clone_reduce(all_expected, expected_result, CLONE_SUM)
    call clone_barrier()
    if (all_expected /= all_sum) then
       if (clone_myid() == 0) then
          write(*,*) clone_myid(), "mothers:", all_sum, all_expected, all_sum - all_expected
       end if
       call printit("face_check mothers", (all_sum ==  all_expected), my_dt)
    end if
    deallocate(values)
  end subroutine faces_check

  subroutine deriv_test(fm, n_iter)
    use clone_lib_module, only: clone_get
    use fakemesh, only: fakemesh_t
    use gradient_types, only: kode_vel
    use matdefcm, only: nummat
    use gradient_types, only: gradient_prop_t
    use interface_types, only: interface_option_t
    use my_derivatives, only: derivatives_common_split
#ifdef ENABLE_VTUNE  
    use ittnotify
#endif
    implicit none
    type(fakemesh_t) :: fm
    integer, intent(in) :: n_iter
    integer :: numitr
    type(gradient_prop_t) :: gradp
    type(interface_option_t) :: intopt
    logical, allocatable :: noslope_cell(:) ! Left unallocated 
    integer :: iMat, iIter
    real(REAL64) :: my_dt
    real(REAL64), allocatable, dimension(:,:)   :: value_cloned
    integer :: iVar
    
    my_dt = now()
#ifdef ENABLE_VTUNE  
    call itt_resume()
#endif

    ASSOCIATE ( &
         m => fm%m, &
         cells => fm%m%cells, &
         core => fm%core, &
         numcell => fm%m%cells%numcell &
         )

      numitr = 1
      gradp%numrho = 1
      gradp%numrho_fvol = 1
      allocate(value_cloned(m%cells%numcell_clone, m%sim%numvel))
      do iVar = 1, m%sim%numvel
         call clone_get(value_cloned(:, iVar))
      end do
      
      do iIter = 1, n_iter
         call derivatives_common_split( &
              m%sim, m, &
              fm%frac_core, fm%core, &
              gradp, intopt, &
              cells%numcell_clone, gradp%numrho, m%sim%numvel,   &
              kode_vel, noslope_cell,  &
              core%deriv_velocity(1:cells%numcell_clone,1:m%sim%numdim,1:m%sim%numvel), &
              .true., &
              invalue = core%cell_velocity(1:cells%numcell_clone,1:m%sim%numvel),&
              value_cloned = value_cloned &
              )
      end do
      deallocate(value_cloned)
    END ASSOCIATE
#ifdef ENABLE_VTUNE  
    call itt_pause()
#endif
    my_dt = now() - my_dt
    call printit("deriv_test", .true., my_dt)

  end subroutine deriv_test
    
end module tests
    
