
program test
  use fakemesh
  use mesh_state_types
  type(fakemesh_t) :: fm
  type(mesh_state_frac_core_t) :: frac_core
  character(len=4096) :: fname
  character(len=4096) :: arg
  integer :: nprocs

  call GET_COMMAND_ARGUMENT(1, fname)
  call GET_COMMAND_ARGUMENT(2, arg)
  read(arg,*) nprocs
  call fm%init_from_PIO(trim(fname))
  write(*,*) 'releasing'
  call fm%release_PIO()
  
  write(*,*) 'woohoo', trim(fname)
end program test
