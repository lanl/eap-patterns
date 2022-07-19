
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
  use fakemesh
  use mesh_state_types
  use clone_lib_module, only: clone_exit, clone_myid, clone_nprocs
  implicit none
  type(fakemesh_t) :: fm
  type(mesh_state_frac_core_t) :: frac_core
  character(len=4096) :: fname
  character(len=4096) :: arg
  integer :: nprocs, myid

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
  if (myid == 0) write(*,*) 'releasing'
  call fm%release_PIO()
  call clone_exit()
  
  if (myid == 0) write(*,*) 'woohoo', trim(fname)
end program test
