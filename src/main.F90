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
  use fakemesh
  use binreader, only: binfile_verify_signature
  use mesh_state_types
  use clone_lib_module, only: clone_exit, clone_myid, clone_nprocs, &
       clone_reduce, CLONE_SUM, clone_barrier
  use tests, only: test_driver
#ifdef ENABLE_VTUNE  
  use ittnotify
#endif
  implicit none
  type(fakemesh_t) :: fm
  type(mesh_state_frac_core_t) :: frac_core
  character(len=4096) :: fname
  character(len=4096) :: arg
  integer :: nprocs, myid
  integer :: n_iter
  real(REAL64) :: my_result, expected_result
  integer(INT64) :: total_numtop, local_numtop
  real(REAL64) :: t0, dt

#ifdef ENABLE_VTUNE  
  call itt_pause()
#endif
  
  ! Get the filename
  call GET_COMMAND_ARGUMENT(1, fname)

  ! Quick check on the file type
  if (.not. binfile_verify_signature(fname)) then
     write(*,*) '_______ERROR: ', trim(fname), ' is not an EAP-bin file'
  else
#ifdef ENABLE_MPI
     call fm%init(trim(fname))
     myid = clone_myid()
     nprocs = clone_nprocs()
#else
     call GET_COMMAND_ARGUMENT(2, arg)
     if ( len_trim(arg) > 0) then
        read(arg,*) nprocs
     else
        nprocs = 1
     end if
     myid = 0
     call fm%init(trim(fname), nprocs, myid)
#endif

     n_iter = 1
     call test_driver(fm, n_iter)

     call clone_barrier()
     if (myid == 0) write(*,*) 'releasing'
     call fm%release()
     call clone_exit()
  endif
end program test
