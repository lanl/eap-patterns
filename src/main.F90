
subroutine testme(fm)
  use sim_types, only : sim_info_t
  use mesh_types, only : mesh_t
  use mesh_state_types, only : mesh_state_frac_core_t, mesh_state_core_t
  use gradient_types,         only : gradient_prop_t
  use interface_types,        only : interface_option_t
  use fakemesh
  implicit none
  type(fakemesh_t) :: fm
  type(mesh_state_frac_core_t) :: frac_core
  
  ! call derivatives_common_splits(fm%m%sim, fm%m, &
  !        frac_core, core, &
  !        gradp, intopt, &
  !        cell_dim, numitr, nvec, kode, &
  !        noslope_cell, deriv, do_fincom, faceval, &
  !        deriv_weight, invalue, value_cloned, do_special)
  !     class(sim_info_t), intent(in) :: sim
  !     type(mesh_t), intent(in) :: mesh
  !     type(mesh_state_frac_core_t), intent(in) :: frac_core
  !     type(mesh_state_core_t), intent(in) :: core
  !     type(gradient_prop_t), intent(in) :: gradp
  !     type(interface_option_t), intent(in) :: intopt
  !     integer,     intent(in) :: cell_dim
  !     integer,     intent(in) :: numitr
  !     integer,     intent(in) :: nvec
  !     integer,     intent(in) :: kode(:,:)
  !     logical,     intent(in), allocatable :: noslope_cell(:)
  !     logical,     intent(in) :: do_fincom
  !     logical,     intent(in), optional :: do_special
end subroutine testme
program test
  use iso_fortran_env, only: REAL64, INT64
  use pio_interface, only: pio_now
  use fakemesh
  use mesh_state_types
  use clone_lib_module, only: clone_exit, clone_myid, clone_nprocs, clone_reduce, CLONE_SUM
  use test_cells, only: test_sum
  implicit none
  type(fakemesh_t) :: fm
  type(mesh_state_frac_core_t) :: frac_core
  character(len=4096) :: fname
  character(len=4096) :: arg
  integer :: nprocs, myid
  integer :: n_iter
  real(REAL64) :: my_result, expected_result
  real(REAL64), allocatable :: values(:)
  integer(INT64) :: total_numtop, local_numtop
  real(REAL64) :: t0, dt

  ! Get the filename
  call GET_COMMAND_ARGUMENT(1, fname)
  
#ifdef EP_MPI

  call fm%init_from_PIO(trim(fname))
  myid = clone_myid()
  nprocs = clone_nprocs()
#else
  call GET_COMMAND_ARGUMENT(2, arg)
  read(arg,*) nprocs
  call fm%init_from_PIO(trim(fname), nprocs, myid)
#endif
  ASSOCIATE(m => fm%m)
    ! Calculate the expected results
    ! Simple_test returns the n_iter * total_numtop
    n_iter = 5
    local_numtop = m%levels%numtop
    call clone_reduce(total_numtop, local_numtop, CLONE_SUM)
    expected_result = real(n_iter,kind=REAL64) * real(total_numtop,kind=REAL64)

    allocate(values(m%cells%numcell))
    values = 1

    t0 = pio_now()
    my_result = test_sum(m, values, n_iter)
    dt = pio_now() - t0
    if (myid == 0) then
       if (my_result /= expected_result) then
          write(*,*) 'Wrong result, expected:', expected_result, ' but got ', my_result, 'in t=', dt
       else
          write(*,*) 'simple test time=', dt
       end if
    end if

  if (myid == 0) write(*,*) 'releasing'
  call fm%release_PIO()
  call clone_exit()



  END ASSOCIATE
end program
