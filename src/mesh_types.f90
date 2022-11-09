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

module timer_module
  implicit none
  public
contains
  subroutine timerset(val, key)
    ! Insert your timer calls here
    implicit none
    logical, intent(in) :: val
    character*(*), intent(in) :: key
    return
  end subroutine timerset
end module timer_module


module mesh_state_types
  use iso_c_binding
  use iso_fortran_env, only : REAL64

  implicit none

  public

  type mesh_state_core_t
     real(REAL64), contiguous, pointer :: rho(:) => null()
     real(REAL64), contiguous, pointer :: cell_velocity(:,:) => null()
     real(REAL64), contiguous, pointer :: deriv_velocity(:,:,:) => null()
  end type mesh_state_core_t
  
  type mesh_state_frac_var_t
     integer :: nmat
     integer :: ncells
     real(REAL64), pointer, dimension(:,:) :: obj
  end type mesh_state_frac_var_t
  
  type mesh_state_frac_core_t
     type(mesh_state_frac_var_t) :: mass
     type(mesh_state_frac_var_t) :: vol       !  Needed
     type(mesh_state_frac_var_t) :: eng

     ! procedure(frac_core_io_int), pointer, nopass :: io => null()
     ! procedure(frac_core_movedn_int), pointer, nopass :: movedn => null()
     ! procedure(frac_core_clip_tiny_volumes_int), pointer, nopass :: clip_tiny_volumes => null()
     ! procedure(frac_core_deallocate_int), pointer, nopass :: deallocate => null()
     ! procedure(frac_core_is_allocated_int), pointer, nopass :: is_allocated => null()
  end type mesh_state_frac_core_t

end module mesh_state_types

module matdefcm
  public
  integer :: nummat = 1
end module matdefcm

module fixed_values_module
  use iso_fortran_env, only : REAL64
  implicit none
  public
  real(REAL64), parameter :: minimum_fraction = 1.00021e-12_REAL64

end module fixed_values_module

module interface_types
  public
  type interface_option_t
     integer :: interface_option = 0
     integer :: vof_multimat_treatment = 1
     logical :: flatten_interface_vel = .false.
     logical, pointer :: is_vof_mat(:) => null()
  end type interface_option_t
end module interface_types

module sim_types
  public
  type ::  sim_info_t
     integer :: numdim=3
     integer :: numvel=3
  end type sim_info_t
end module sim_types

module var_wrapper_class
  use iso_fortran_env, only : REAL64
  public
  type :: var_wrapper
  end type var_wrapper
contains
  subroutine vw_set(vw, deriv, ncell, ndim, nvec)
    type(var_wrapper), intent(inout) :: vw
    real(REAL64), dimension(:,:,:), intent(in) :: deriv
    integer, intent(in) :: ncell
    integer, intent(in) :: ndim
    integer, intent(in) :: nvec
  end subroutine vw_set
end module var_wrapper_class

module gradient_types
  public
  integer, parameter :: KODE_LEN = 3
  integer, parameter :: kode_zero(KODE_LEN) = 0
  integer, parameter :: kode_vel(KODE_LEN,KODE_LEN)  = &
       reshape ( &
       [ -1, -10, -10, & 
       -10, -1,  -10, &
       -10, -10, -1 ], [KODE_LEN,KODE_LEN])
  type gradient_prop_t
     integer :: shock_detector
     integer :: numrho
     integer :: numrho_fvol
  end type gradient_prop_t
end module gradient_types


module mesh_state_cell_accessors
  public
contains
  pure subroutine CV_zero_for_cells_ary(weight,mixed_count,mixed_list)
    use iso_fortran_env, only : REAL64
    implicit none
    real(REAL64), intent(in) :: weight(:)
    integer, intent(in) :: mixed_count
    integer, intent(inout) :: mixed_list(:)
  end subroutine CV_zero_for_cells_ary
  pure function CV_get_for_cells_ary(var,cell_count,cell_list)
    use iso_fortran_env, only : REAL64
    implicit none
    real(REAL64), dimension(:), pointer, intent(in) :: var
    integer,intent(in) :: cell_count, cell_list(:)
    real(REAL64), dimension(:) :: CV_get_for_cells_ary(cell_count)
    CV_get_for_cells_ary =  var(cell_list)
  end function CV_get_for_cells_ary
end module mesh_state_cell_accessors



module mesh_scratch_gravity
  use iso_fortran_env, only : REAL64
  public
  ! ------------------------------------------------------------------------------
  type mesh_scratch_gravity_t
     real(REAL64), dimension(:,:), allocatable :: grav_accel
  end type mesh_scratch_gravity_t
  ! ------------------------------------------------------------------------------
end module mesh_scratch_gravity
module mesh_scratch_gravity_module
  use mesh_scratch_gravity
  public
  type(mesh_scratch_gravity_t) ::  grav_scr
end module mesh_scratch_gravity_module

module mesh_types
  use, intrinsic :: iso_fortran_env, only : INT64, REAL64
  use iso_c_binding
  use sim_types, only: sim_info_t
  use mesh_state_types
  public
  type levels_t
     integer, pointer :: numtop => null()
     integer, pointer :: allnumtop => null()
     integer, dimension(:), pointer :: ltop => null()
     integer, dimension(:), pointer :: alltop => null()
     integer, dimension(:), pointer :: ltop_nv => null()
     integer, dimension(:), pointer :: cell_level => null()
     integer(INT64), dimension(:), pointer :: cell_daughter => null()
  end type levels_t
  type faces_t
     integer, dimension(:), pointer :: face_num => null()
     integer, dimension(:,:), pointer :: face_hi => null()
     integer, dimension(:,:), pointer :: face_lo => null()
     integer, dimension(:,:), pointer :: face_flag => null()
     integer, dimension(:,:), pointer :: face_id => null()
     integer, dimension(:,:,:), pointer :: face_local => null()
  end type faces_t
  type neighbors_t
  end type neighbors_t
  type amr_vars_t
  end type amr_vars_t
  type user_refine_vars_t
  end type user_refine_vars_t
  type cells_t
     integer, pointer :: numcell => null()                          ! needed
     integer(INT64), pointer ::   sum_numcell  => null()
     integer, pointer :: max_numcell => null()
     integer, pointer :: numcell_clone => null()                    ! needed
     integer, pointer :: mxcell => null() 

     integer(INT64), dimension(:), pointer :: cell_address => null()

     logical(c_bool), dimension(:), pointer :: cell_active => null()

     ! set in kidmom_module in check_cell_center
     real(REAL64), dimension(:,:), pointer :: cell_center => null()

     ! Geometric centroid of cell, same for rectangular geometry, differs
     ! for cylindrical and spherical geometry.
     real(REAL64), dimension(:,:), pointer :: cell_position => null()

     real(REAL64), dimension(:,:), pointer :: cell_half => null()

     ! Radial distance from the cell lower/high edge to the centroid.
     real(REAL64), dimension(:,:), pointer :: cell_half_lo => null(), &    ! needed 
          cell_half_hi => null()    ! needed 


     ! cell volume
     real(REAL64), dimension(:), pointer :: vcell                ! needed

     integer, dimension(:), pointer :: global_numcell
     integer(INT64), dimension(:), pointer :: global_base, global_base_old

  end type cells_t
  type mesh_t
     class(cells_t), allocatable :: cells
     class(faces_t), allocatable :: faces
     class(levels_t), allocatable :: levels
     class(neighbors_t), allocatable :: neighbors
     class(sim_info_t), pointer :: sim
     class(amr_vars_t), allocatable :: amr_vars
     type(user_refine_vars_t) :: usref_vars
  end type mesh_t
  contains
    subroutine release_mesh(m)
      use mem_release, only: release
      type(mesh_t) :: m
      if (allocated(m%cells)) then
         ! Release cells
         
         call release(m%cells%numcell)
         call release(m%cells%sum_numcell)
         call release(m%cells%max_numcell)
         call release(m%cells%numcell_clone)
         call release(m%cells%mxcell)
         call release(m%cells%cell_address)
         call release(m%cells%cell_active)
         call release(m%cells%cell_center)
         call release(m%cells%cell_position)
         call release(m%cells%cell_half)
         call release(m%cells%cell_half_lo)
         call release(m%cells%cell_half_hi)
         call release(m%cells%vcell)
         call release(m%cells%global_numcell)
         call release(m%cells%global_base)
         call release(m%cells%global_base_old)
         deallocate(m%cells)
      end if
      if (allocated(m%levels)) then
         ! Release levels
         call release(m%levels%cell_daughter)
         call release(m%levels%cell_level)
         call release(m%levels%numtop)
         call release(m%levels%allnumtop)
         call release(m%levels%ltop)
         call release(m%levels%alltop)
         deallocate(m%levels)
      end if
      if (allocated(m%faces)) then
         ! Release faces
         call release(m%faces%face_num)
         call release(m%faces%face_hi)
         call release(m%faces%face_lo)
         call release(m%faces%face_flag)
         call release(m%faces%face_id)
         call release(m%faces%face_local)
         deallocate(m%faces)
      end if
      if(allocated(m%neighbors)) then
         deallocate(m%neighbors)
      end if
      if(associated(m%sim)) then
         deallocate(m%sim)
      end if
      if(allocated(m%amr_vars)) then
         deallocate(m%amr_vars)
      end if
    end subroutine release_mesh
    subroutine nullify_mesh(m)
      type(mesh_t) :: m
      ! nullify cell members
      nullify(m%cells%numcell)
      nullify(m%cells%sum_numcell)
      nullify(m%cells%max_numcell)
      nullify(m%cells%numcell_clone)
      nullify(m%cells%mxcell)
      nullify(m%cells%cell_address)
      nullify(m%cells%cell_active)
      nullify(m%cells%cell_center)
      nullify(m%cells%cell_position)
      nullify(m%cells%cell_half)
      nullify(m%cells%cell_half_lo)
      nullify(m%cells%cell_half_hi)
      nullify(m%cells%vcell)
      nullify(m%cells%global_numcell)
      nullify(m%cells%global_base)
      nullify(m%cells%global_base_old)

      ! nullify face members
      nullify(m%faces%face_num)
      nullify(m%faces%face_hi)
      nullify(m%faces%face_lo)
      nullify(m%faces%face_flag)
      nullify(m%faces%face_id)
      nullify(m%faces%face_local)
    end subroutine nullify_mesh
end module mesh_types

module mesh_state_accessors
  use iso_fortran_env, only : REAL64
  implicit none
  public
contains
  pure integer function FV_count_mixed_top_cells_range(vol,mesh,nstart,nend,nummat)
    use mesh_state_types, only: mesh_state_frac_var_t
    use mesh_types, only: mesh_t
    implicit none
    type(mesh_state_frac_var_t), intent(in) :: vol
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: nstart, nend, nummat
  end function FV_count_mixed_top_cells_range
  pure function FV_get_for_cells_ary_materials(vol,cell_count,cell_list,vmatnum,vmats)
    use mesh_state_types, only: mesh_state_frac_var_t
    use mesh_types, only: mesh_t
    implicit none
    type(mesh_state_frac_var_t), intent(in) :: vol
    integer,intent(in) :: cell_count, cell_list(:)
    integer,intent(in) :: vmatnum, vmats(:)
    real(REAL64), dimension(:) :: FV_get_for_cells_ary_materials(cell_count, vmatnum)
  end function FV_get_for_cells_ary_materials
end module mesh_state_accessors
