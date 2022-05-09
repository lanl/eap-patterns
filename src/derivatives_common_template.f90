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

  ! Reconstituted to use simpler data structures and fed from
  ! captured data
#if defined(__INTEL_COMPILER) || defined(_SEQUOIA) || defined(CRAY)
#define NOT_GCC
#endif
  
#if defined(_SEQUOIA) && defined(__GNUC__)
#undef NOT_GCC
#endif
  
#ifndef NOT_GCC
  /* The fortran preprocessor is in traditional mode, so token concatentation is like below:*/
#define stringify(s1) "s1"
#else
#define stringify(s1) #s1
#endif

#ifdef CRAY
#define stringify(s1) "s1"
#endif

#define stringjoin(s1,s2) (stringify(s1) // "_" // stringify(s2))
#define TIMER_NAME(key) stringjoin(derivatives_common,key)
#define TIMERSET(val,key)

  ! Define this as needed to get timer
#define TIMERSET_UNCOND(val,key) call timerset(val,TIMER_NAME(key))

#define get_value(l,nm) merge(value_cloned(l,nm),invalue(l,nm),l .gt. mesh%cells%numcell)

  module my_derivatives
    use iso_fortran_env, only : REAL64
    use define_kind, only: ZERO, HALF, ONE, TWO, HI_SIDE, LO_SIDE, mype, iope
    use sim_types, only : sim_info_t
    use mesh_types, only : mesh_t, cells_t, faces_t
    use mesh_state_types, only : mesh_state_frac_core_t, mesh_state_core_t
    implicit none
    public
    integer, parameter :: NO_DERIV = 0
    integer, parameter :: INT_MM = 1
    integer, parameter :: INT_EMM = 2
    integer, parameter :: INT_LV = 3
    integer, parameter :: INT_NL = 4

    integer :: itrmax = 1
    integer :: method = 0
    integer :: nitr = 1
    logical :: limit_slope = .false.
  contains
      pure function mixed_volume_fraction(fv)
        use fixed_values_module, only : minimum_fraction
        real(REAL64), parameter         :: cvl = minimum_fraction !Threshold below which a material vol. frac. is considered zero
        real(REAL64), parameter         :: cvh = ONE - cvl        !Threshold above which a material vol. frac. is considered one
        real(REAL64), intent(in) :: fv
        logical                  :: mixed_volume_fraction

        mixed_volume_fraction = (fv.gt.cvl .and. fv.lt.cvh)

      end function mixed_volume_fraction
! ------------------------------------------------------------------------------
    subroutine modify_weight_tempo(mesh, frac_core, intopt, weight)
      use mesh_state_accessors, only : fv_count_mixed_top_cells_range
      use mesh_state_cell_accessors, only : cv_zero_for_cells_ary
      use matdefcm,        only : nummat
      use interface_types, only : interface_option_t
      use util          ,  only : global_error
      use timer_module, only : timerset

      ! ----- purpose

      ! ----- flag cells that are vof_cells. Create a compressed linked list of just
      ! ----- cells that will use vof reconstruction.

      ! ----- calling arguments
      type(mesh_t), intent(in) :: mesh
      type(mesh_state_frac_core_t), intent(in) :: frac_core
      type(interface_option_t), intent(in) :: intopt
      real(REAL64), intent(out)    :: weight(:)

      ! ----- local variables

      integer       :: nt, nm, l
      integer       :: thr
      integer       :: nstart, nend, nlen
      integer :: vmats(nummat), vmatnum, vmatpos
      integer, parameter :: step_size = 512
      integer :: counts(step_size)
      integer :: cell_list(step_size), cell_count
      integer :: mixed_list(step_size), mixed_count

      associate ( cells => mesh%cells, &
           levs => mesh%levels)

        TIMERSET(.true., modify_weight_tempo)

        do l  = 1,mesh%cells%numcell_clone
           weight(l) = ONE
        enddo

        if(intopt%interface_option.ge.3) then
           vmatnum = 0
           do nm = 1,nummat
              if(intopt%is_vof_mat(nm)) then
                 vmatnum = vmatnum + 1
                 vmats(vmatnum) = nm
              endif
           enddo

           do nstart = 1,mesh%levels%numtop,step_size
              nend = min(nstart+step_size-1,mesh%levels%numtop)
              nlen = nend-nstart+1
              counts(1:nlen) = FV_count_mixed_top_cells_range(frac_core%vol,mesh,nstart,nend,nummat)

              cell_count = 0
              if (intopt%vof_multimat_treatment.eq.1) then
                 do nt = 1,nlen
                    if (counts(nt).ge.2) then
                       cell_count = cell_count + 1
                       cell_list(cell_count) = levs%ltop_nv(nstart+nt-1)
                    endif
                 enddo
              elseif (intopt%vof_multimat_treatment.eq.2) then
                 do nt = 1,nlen
                    if (counts(nt).eq.2) then
                       cell_count = cell_count + 1
                       cell_list(cell_count) = levs%ltop_nv(nstart+nt-1)
                    endif
                 enddo
              endif

              call collect_vof_mat_mixed(cells,frac_core,cell_count,cell_list,vmatnum,vmats,mixed_count,mixed_list)
              call CV_zero_for_cells_ary(weight,mixed_count,mixed_list(1:mixed_count))
           enddo
        endif

        TIMERSET(.false., modify_weight_tempo)

      end associate

    end subroutine modify_weight_tempo
    ! ------------------------------------------------------------------------------
    subroutine collect_vof_mat_mixed(cells, frac_core, cell_count,cell_list, &
         & vmatnum, vmats, mixed_count, mixed_list)

      use mesh_state_accessors, only : fv_get_for_cells_ary_materials
      use mesh_state_cell_accessors, only : cv_get_for_cells_ary

      class(cells_t),intent(in) :: cells
      type(mesh_state_frac_core_t), intent(in) :: frac_core
      integer,intent(in) :: cell_count, cell_list(:)
      integer,intent(in) :: vmatnum, vmats(:)
      integer, intent(out) :: mixed_count, mixed_list(:)

      real(REAL64) :: vol(cell_count,vmatnum), cvols(cell_count)
      integer :: nt, l
      integer :: vmatpos
      logical :: is_mixed
      real(REAL64)  :: relvol, copy_fvol

      vol(1:cell_count,1:vmatnum) = &
           & FV_get_for_cells_ary_materials(frac_core%vol,cell_count,cell_list,vmatnum,vmats)
      cvols(1:cell_count) = CV_get_for_cells_ary(cells%vcell,cell_count,cell_list)
      mixed_count = 0

      do nt = 1,cell_count
         l = cell_list(nt)

         is_mixed = .false.
         do vmatpos = 1,vmatnum
            relvol = vol(nt,vmatpos) / cvols(nt)
            copy_fvol  = max(ZERO,min(ONE,relvol))
            if (mixed_volume_fraction(copy_fvol)) then  ! and vofmat...
               is_mixed = .true.
               exit
            endif
         enddo
         if (is_mixed) then
            mixed_count = mixed_count + 1
            mixed_list(mixed_count) = l
         endif
      enddo
    end subroutine collect_vof_mat_mixed
    subroutine inside_com3b(sim, mesh, dir, nm, cell_val_mnmx_hilo, kode, &
         & cell_value_hilo, do_special, do_pressure, core, faceval )

      ! non-optional scalars and derived types
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      integer,      intent(in)     :: dir, nm

      ! non-optional arrays
      real(REAL64), intent(inout), dimension(:,:) :: cell_val_mnmx_hilo !(mesh%cells%numcell_clone,4)
      integer,      intent(in)     :: kode(:,:)
      real(REAL64), intent(inout) :: cell_value_hilo(:,:,:)

      ! optional scalars and derived types
      logical,      intent(in), optional       :: do_special
      logical,      intent(in), optional       :: do_pressure
      type(mesh_state_core_t), optional, intent(in) :: core

      ! optional arrays
      !       for inflow bndy cond'ns
      real(REAL64), intent(in), optional     :: faceval(:,:)

      ! ----- local variables

      integer :: loop, n, lhi, llo, lmi, lmo
      real(REAL64) :: face_value

      associate (cells => mesh%cells, &
           faces => mesh%faces)

        TIMERSET(.true., inside_com3b)

        ! .....    calculate the area weighted average face values
        do loop = 1,faces%face_num(dir)

           if (faces%face_id(loop,dir) .gt. 2) then

              if (present(do_special)) then
                 ! Do the HI_SIDE part
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    lhi = faces%face_local(n,HI_SIDE,dir)
                    llo = faces%face_local(n,LO_SIDE,dir)
                    face_value = ZERO

                    if (core%rho(llo).gt.ZERO .or. core%rho(lhi).gt.ZERO) then
                       if (do_pressure .and. cell_value_hilo(llo,1,nm)       &
                            & * cell_value_hilo(lhi,2,nm).le.ZERO) then
                          ! .....  this coding addresses the hot spot problem
                          face_value = (cells%cell_half_lo(lhi,dir)*core%rho(lhi) *  &
                               & cell_value_hilo(llo,2,nm)         &
                               & +  cells%cell_half_hi(llo,dir)*core%rho(llo) *  &
                               & cell_value_hilo(lhi,1,nm))        &
                               & / (cells%cell_half_lo(lhi,dir)*core%rho(lhi)    &
                               & +  cells%cell_half_hi(llo,dir)*core%rho(llo))
                       else
                          face_value = (cells%cell_half_lo(lhi,dir) *           &
                               & cell_value_hilo(llo,2,nm)         &
                               & +  cells%cell_half_hi(llo,dir) *           &
                               & cell_value_hilo(lhi,1,nm))        &
                               & / (cells%cell_half_lo(lhi,dir)             &
                               & +  cells%cell_half_hi(llo,dir))
                       endif ! do_pressure ...
                    endif ! core%rho
                    cell_val_mnmx_hilo(lhi,1) =                          &
                         & min(cell_val_mnmx_hilo(lhi,1), face_value)
                    cell_val_mnmx_hilo(lhi,2) =                          &
                         & max(cell_val_mnmx_hilo(lhi,2), face_value)
                 enddo ! n

                 ! do the LO_SIDE part
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    llo = faces%face_local(n,LO_SIDE,dir)
                    lhi = faces%face_local(n,HI_SIDE,dir)
                    face_value = ZERO

                    if (core%rho(llo).gt.ZERO .or. core%rho(lhi).gt.ZERO) then
                       if (do_pressure .and. cell_value_hilo(llo,1,nm)       &
                            & *cell_value_hilo(lhi,2,nm).le.ZERO) then
                          ! .....  this coding addresses the hot spot problem
                          face_value = (cells%cell_half_lo(lhi,dir)*core%rho(lhi) *  &
                               & cell_value_hilo(llo,2,nm)         &
                               & +  cells%cell_half_hi(llo,dir)*core%rho(llo) *  &
                               & cell_value_hilo(lhi,1,nm))        &
                               & / (cells%cell_half_lo(lhi,dir)*core%rho(lhi)    &
                               & +  cells%cell_half_hi(llo,dir)*core%rho(llo))
                       else
                          face_value = (cells%cell_half_lo(lhi,dir) *           &
                               & cell_value_hilo(llo,2,nm)         &
                               & +  cells%cell_half_hi(llo,dir) *           &
                               & cell_value_hilo(lhi,1,nm))        &
                               & / (cells%cell_half_lo(lhi,dir)             &
                               & +  cells%cell_half_hi(llo,dir))
                       endif
                    endif
                    cell_val_mnmx_hilo(llo,3) =                          &
                         & min(cell_val_mnmx_hilo(llo,3), face_value)
                    cell_val_mnmx_hilo(llo,4) =                          &
                         & max(cell_val_mnmx_hilo(llo,4), face_value)
                 enddo ! n

              else !not present(do_special)
                 !       Do the HI_SIDE first
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    lhi = faces%face_local(n,HI_SIDE,dir)
                    llo = faces%face_local(n,LO_SIDE,dir)
                    face_value = (cells%cell_half_lo(lhi,dir)*cell_value_hilo(llo,2,nm)  &
                         & +cells%cell_half_hi(llo,dir)*cell_value_hilo(lhi,1,nm)) &
                         & /(cells%cell_half_lo(lhi,dir)+cells%cell_half_hi(llo,dir))
                    cell_val_mnmx_hilo(lhi,1) = min(cell_val_mnmx_hilo(lhi,1), face_value)
                    cell_val_mnmx_hilo(lhi,2) = max(cell_val_mnmx_hilo(lhi,2), face_value)
                 enddo ! n

                 !       Do the LO_SIDE
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    llo = faces%face_local(n,LO_SIDE,dir)
                    lhi = faces%face_local(n,HI_SIDE,dir)
                    face_value = (cells%cell_half_lo(lhi,dir)*cell_value_hilo(llo,2,nm)  &
                         & +cells%cell_half_hi(llo,dir)*cell_value_hilo(lhi,1,nm)) &
                         & /(cells%cell_half_lo(lhi,dir)+cells%cell_half_hi(llo,dir))
                    cell_val_mnmx_hilo(llo,3) = min(cell_val_mnmx_hilo(llo,3), face_value)
                    cell_val_mnmx_hilo(llo,4) = max(cell_val_mnmx_hilo(llo,4), face_value)
                 enddo ! n
              endif ! present(do_special)

           else if (faces%face_id(loop,dir) .eq. 2) then

              if(kode(dir,nm).eq.-2)then
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    lmo = faces%face_local(n,LO_SIDE,dir)
                    face_value = faceval(n,dir)
                    cell_val_mnmx_hilo(lmo,3) = min(cell_val_mnmx_hilo(lmo,3), face_value)
                    cell_val_mnmx_hilo(lmo,4) = max(cell_val_mnmx_hilo(lmo,4), face_value)
                 enddo ! n
              else
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    llo = faces%face_local(n,LO_SIDE,dir)
                    face_value = cell_value_hilo(llo,2,nm)
                    cell_val_mnmx_hilo(llo,3) = min(cell_val_mnmx_hilo(llo,3), face_value)
                    cell_val_mnmx_hilo(llo,4) = max(cell_val_mnmx_hilo(llo,4), face_value)
                 enddo ! n
              endif

           else if (faces%face_id(loop,dir) .eq. 1) then

              if(kode(dir,nm).eq.-2)then
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    lmi = faces%face_local(n,HI_SIDE,dir)
                    face_value = faceval(n,dir)
                    cell_val_mnmx_hilo(lmi,1) = min(cell_val_mnmx_hilo(lmi,1), face_value)
                    cell_val_mnmx_hilo(lmi,2) = max(cell_val_mnmx_hilo(lmi,2), face_value)
                 enddo ! n
              else
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    lhi = faces%face_local(n,HI_SIDE,dir)
                    face_value = cell_value_hilo(lhi,1,nm)
                    cell_val_mnmx_hilo(lhi,1) = min(cell_val_mnmx_hilo(lhi,1), face_value)
                    cell_val_mnmx_hilo(lhi,2) = max(cell_val_mnmx_hilo(lhi,2), face_value)
                 enddo ! n
              endif

           endif ! face_id

        enddo ! loop

        TIMERSET(.false., inside_com3b)

      end associate

    end subroutine inside_com3b
    ! ------------------------------------------------------------------------------
    subroutine inside_com3e(sim, mesh, itr, nm, limit_slope, cell_deriv_hilo, &
         & deriv_mm )

      ! non-optional scalars and derived types
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      integer,      intent(in)     :: itr, nm
      logical, intent(in) :: limit_slope

      ! non-optional arrays
      real(REAL64), intent(inout), dimension(:,:) :: cell_deriv_hilo !(mesh%cells%numcell_clone,2)
      real(REAL64), intent(inout) :: deriv_mm(:,:)

      ! ----- local variables
      integer      :: nt, l

      associate (cells => mesh%cells, &
           & faces => mesh%faces)

        TIMERSET(.true., inside_com3e)

        ! .....       
        if (limit_slope .and. itr.eq.1) then
           do nt = 1,mesh%levels%allnumtop
              l   = mesh%levels%alltop(nt)
              deriv_mm(l,nm) = &
                   &   max(ZERO, min(cell_deriv_hilo(l,1),cell_deriv_hilo(l,2))) &
                   & + min(ZERO, max(cell_deriv_hilo(l,1),cell_deriv_hilo(l,2)))
           enddo ! nt
        endif

        TIMERSET(.false., inside_com3e)

      end associate

    end subroutine inside_com3e
    ! ------------------------------------------------------------------------------
    subroutine inside_com3f(sim, mesh, dir, nm, numitr, method, &
         & lv_weight, deriv, cell_deriv_hilo, cell_deriv_hilo_all_dir )

      ! non-optional scalars and derived types
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      integer, intent(in) :: dir, nm, numitr, method
      real(REAL64), intent(in)     :: lv_weight

      ! non-optional arrays
      real(REAL64), intent(out)    :: deriv(:,:,:)
      real(REAL64), intent(inout), dimension(:,:) :: cell_deriv_hilo !(mesh%cells%numcell_clone,2)
      real(REAL64), intent(inout) :: cell_deriv_hilo_all_dir(:,:,:,:) !JV

      ! ----- local variables
      integer      :: nt, l

      associate (cells => mesh%cells, &
           faces => mesh%faces)

        TIMERSET(.true., inside_com3f)

        ! .....       calculate the cell derivative
        select case (method)
        case (INT_MM)  
           do nt = 1,mesh%levels%allnumtop
              l   = mesh%levels%alltop(nt)
              deriv(l,dir,nm) = max(ZERO,min(cell_deriv_hilo(l,1),cell_deriv_hilo(l,2))) &
                   & +min(ZERO,max(cell_deriv_hilo(l,1),cell_deriv_hilo(l,2)))
           enddo ! nt

        case (INT_EMM)  

           do nt = 1,mesh%levels%allnumtop
              l   = mesh%levels%alltop(nt)
              deriv(l,dir,nm) = median(cell_deriv_hilo(l,1), &
                   & cell_deriv_hilo(l,2), &
                   & -(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))
           enddo ! nt

        case (INT_LV) 

           if( numitr.eq.7 .or. numitr.eq.8 ) then !JV store left and right derivatives
              do nt = 1,mesh%levels%allnumtop
                 l   = mesh%levels%alltop(nt)
                 cell_deriv_hilo_all_dir(dir,l,1:2,nm) = cell_deriv_hilo(l,1:2)
                 deriv(l,dir,nm) = (max(ZERO, min(lv_weight*cell_deriv_hilo(l,1), &
                      & lv_weight*cell_deriv_hilo(l,2), &
                      & HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))) &
                      & +min(ZERO, max(lv_weight*cell_deriv_hilo(l,1), &
                      & lv_weight*cell_deriv_hilo(l,2), &
                      & HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))))
              enddo ! nt
           else 
              do nt = 1,mesh%levels%allnumtop
                 l   = mesh%levels%alltop(nt)
                 deriv(l,dir,nm) = (max(ZERO, min(lv_weight*cell_deriv_hilo(l,1), &
                      & lv_weight*cell_deriv_hilo(l,2), &
                      & HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))) &
                      & +min(ZERO, max(lv_weight*cell_deriv_hilo(l,1), &
                      & lv_weight*cell_deriv_hilo(l,2), &
                      & HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))))
              enddo ! nt
           endif

        case (INT_NL) ! --- turn off limiters

           do nt = 1,mesh%levels%allnumtop
              l   = mesh%levels%alltop(nt)
              deriv(l,dir,nm) = HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2))
           enddo ! nt

        end select

        TIMERSET(.false., inside_com3f)

      end associate

    end subroutine inside_com3f
    ! ------------------------------------------------------------------------------
    subroutine inside_com3g(mesh, dir, numvec, itr, itrmax, limit_slope, &
         & deriv, deriv_mm, gradp, deriv_weight )

      use gradient_types,        only : gradient_prop_t
      use interface_types,       only : interface_option_t

      ! non-optional scalars and derived types
      type(mesh_t), intent(in) :: mesh
      integer,      intent(in)     :: dir, numvec, itr, itrmax
      logical, intent(in)    :: limit_slope

      ! non-optional arrays
      real(REAL64), intent(out)    :: deriv(:,:,:)
      real(REAL64), intent(inout) :: deriv_mm(:,:)


      ! optional scalars and derived types
      type(gradient_prop_t), optional, intent(in) :: gradp

      ! optional arrays
      real(REAL64), intent(in), optional :: deriv_weight(:)


      ! ----- local variables
      integer      :: nm, nt, l

      associate (cells => mesh%cells, &
           faces => mesh%faces)

        TIMERSET(.true., inside_com3g)

        ! ..... limit the slope
        if (itr .eq. itrmax) then
           if (present(gradp)) then
              if (limit_slope .and. present(deriv_weight) .and. gradp%shock_detector.ge.1) then
                 do nm   = 1,numvec
                    do nt = 1,mesh%levels%allnumtop
                       l   = mesh%levels%alltop(nt)
                       deriv(l,dir,nm) = deriv_mm(l,nm) +          &
                            & deriv_weight(l)*(deriv(l,dir,nm)-deriv_mm(l,nm))
                    enddo ! nt
                 enddo ! nm
              endif
           endif
        endif


        TIMERSET(.false., inside_com3g)

      end associate

    end subroutine inside_com3g
    ! ------------------------------------------------------------------------------
    subroutine inside_com3h(sim, mesh, dir, numitr, numvec, &
         & kode, deriv, frac_core, intopt )

      use interface_types,       only : interface_option_t

      ! non-optional scalars and derived types
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      integer,      intent(in)     :: dir, numitr, numvec

      ! non-optional arrays
      integer,      intent(in)     :: kode(:,:)
      real(REAL64), intent(out)    :: deriv(:,:,:)

      ! optional scalars and derived types
      type(mesh_state_frac_core_t), optional, intent(in) :: frac_core
      type(interface_option_t), optional, intent(in) :: intopt

      ! ----- local variables
      integer      :: nm, l, dnr, n
      real(REAL64) :: vof_weight(mesh%cells%numcell_clone)

      associate (cells => mesh%cells, &
           faces => mesh%faces)

        TIMERSET(.true., inside_com3h)

        if (dir .ge. sim%numdim) then
           if (present(intopt)) then
              ! ..... limit the slope to ZERO at interface fronts
              if (intopt%interface_option.ge.3) then

                 ! ..... Get where are the interface cells
                 do l = 1,mesh%cells%numcell_clone
                    vof_weight(l) = ONE
                 enddo
                 if (numitr.ge.0) call modify_weight_tempo(mesh, frac_core, intopt, vof_weight)
                 !-mf mar07
                 do nm    = 1,numvec
                    do dnr = 1,sim%numdim
                       if ((kode(dnr,nm).eq.-1.or.kode(dnr,nm).eq.-10 .or.                  &
                            & kode(dnr,nm).eq.30.or.kode(dnr,nm).eq.2).and. &
                            & .not. intopt%flatten_interface_vel) then
                          !-vel, pres, grav do nothing
                       else
                          do n  = 1,mesh%levels%allnumtop
                             l  = mesh%levels%alltop(n)
                             deriv(l,dnr,nm) = deriv(l,dnr,nm)*vof_weight(l)
                          enddo ! nt
                       endif
                    enddo
                 enddo
              endif
           endif
        endif

        TIMERSET(.false., inside_com3h)

      end associate

    end subroutine inside_com3h
    ! ------------------------------------------------------------------------------
    pure real(REAL64) function median(x, y, z)
      implicit none

      real(REAL64), intent(in)  :: x, y, z

      real(REAL64) :: s

      s = HALF*(sign(ONE, (y-x))+sign(ONE, (z-x)))
      median = x + s*min(abs(y-x), abs(z-x))
    end function median

    ! ------------------------------------------------------------------------------
    subroutine inside_com1_split(sim, mesh, core, numvec,&
         kode, deriv, cell_val_flcl, invalue, value_cloned)
      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t
      use mesh_state_types,    only : mesh_state_core_t
      use mesh_scratch_gravity_module, only : grav_scr
      use util,                only : global_error
      use timer_module, only : timerset

      ! ----- calling arguments
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      type(mesh_state_core_t), intent(in) :: core
      integer,      intent(in)    :: numvec

      integer,      intent(in)    :: kode(:,:)

      real(REAL64), intent(out)   :: deriv(:,:,:)
      real(REAL64), intent(inout), dimension(:,:,:) :: cell_val_flcl
      real(REAL64), intent(in), dimension(:,:) :: invalue
      real(REAL64), intent(in), optional, dimension(:,:) :: value_cloned

      ! ----- local variables

      integer      :: nm, dir, nt, l, loop, n, llo, lhi

      ! ----- code

      associate (cells => mesh%cells, &
           faces => mesh%faces)

        TIMERSET(.true., inside_com1)

        ! .... initialize derivatives
        !      The variable rho must have been communicated if kode(dir,nm).eq.2
        !      and method.ne.NO_DERIV (see do lopp below)

        ! A faster initialization !
        ! deriv = ZERO
        ! do nm = 1,numvec
        !    do dir = 1,sim%numdim
        !       if (method.ne.NO_DERIV .and. kode(dir,nm).eq.2) then
        !          deriv(1:mesh%cells%numcell_clone,dir,nm) = core%rho(1:mesh%cells%numcell_clone) * grav_scr%grav_accel(1:mesh%cells%numcell_clone,dir)
        !       endif
        !    enddo ! dir
        ! enddo ! nm

        do nm = 1,numvec
           do dir = 1,sim%numdim
              if (method.ne.NO_DERIV .and. kode(dir,nm).eq.2) then
                 do l =  1,mesh%cells%numcell_clone
                    deriv(l,dir,nm) = core%rho(l) * grav_scr%grav_accel(l,dir)
                 enddo
              else
                 deriv(1:mesh%cells%numcell_clone,dir,nm) = ZERO
              endif
           enddo ! dir
        enddo ! nm

        ! ..... calculate non-zero derivatives

        if (itrmax .ge. 1) then
           cell_val_flcl(1:mesh%cells%numcell_clone,1:2,1:numvec) = ZERO

           do nm = 1,numvec

              if (mesh%levels%allnumtop.gt.0) then
                 do n = 1,mesh%levels%allnumtop
                    l = mesh%levels%alltop(n)
                    cell_val_flcl(l,1,nm) = get_value(l,nm)
                    cell_val_flcl(l,2,nm) = get_value(l,nm)
                 enddo

                 do dir = 1,sim%numdim
                    select case (kode(dir,nm))
                    case (-1, -10, -20)

                       do loop = 1,faces%face_num(dir)

                          if (kode(dir,nm).eq.-1 .or. &
                               (kode(dir,nm).eq.-10 .and. faces%face_flag(loop,dir).eq.2) .or. &
                               (kode(dir,nm).eq.-20 .and. faces%face_flag(loop,dir).eq.0)) then

                             if (faces%face_id(loop,dir).eq.1) then
                                do n = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                                   lhi = faces%face_local(n,HI_SIDE,dir)
                                   cell_val_flcl(lhi,1,nm) = min(cell_val_flcl(lhi,1,nm), ZERO)
                                   cell_val_flcl(lhi,2,nm) = max(cell_val_flcl(lhi,2,nm), ZERO)
                                enddo ! n
                             else if (faces%face_id(loop,dir).eq.2) then
                                do n = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                                   llo = faces%face_local(n,LO_SIDE,dir)
                                   cell_val_flcl(llo,1,nm) = min(cell_val_flcl(llo,1,nm), ZERO)
                                   cell_val_flcl(llo,2,nm) = max(cell_val_flcl(llo,2,nm), ZERO)
                                enddo ! n
                             endif ! face_id

                          endif
                       enddo ! loop

                    case (2)
                       if (method.ne.NO_DERIV) then
                          do nt = 1,mesh%levels%allnumtop
                             l   = mesh%levels%alltop(nt)
                             if(deriv(l,dir,nm).gt.ZERO)then
                                cell_val_flcl(l,1,nm) = cell_val_flcl(l,1,nm) &
                                     & - deriv(l,dir,nm)*cells%cell_half_lo(l,dir)
                                cell_val_flcl(l,2,nm) = cell_val_flcl(l,2,nm) &
                                     & + deriv(l,dir,nm)*cells%cell_half_hi(l,dir)
                             else if(deriv(l,dir,nm).lt.ZERO)then
                                cell_val_flcl(l,2,nm) = cell_val_flcl(l,2,nm) &
                                     & - deriv(l,dir,nm)*cells%cell_half_lo(l,dir)
                                cell_val_flcl(l,1,nm) = cell_val_flcl(l,1,nm) &
                                     & + deriv(l,dir,nm)*cells%cell_half_hi(l,dir)
                             endif
                          enddo ! nt
                       endif
                    end select
                 enddo ! dir
              endif
           enddo ! nm

           !       Not necessary to communicate ceiling and floor here because they
           !       have been estimated in the ghost cells.

           do nm = 1,numvec
              if (mesh%levels%allnumtop.gt.0) then
                 do dir = 1,sim%numdim
                    do loop = 1,faces%face_num(dir)

                       if (faces%face_id(loop,dir).gt.2) then

                          !                   HI_SIDE
                          do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                             lhi = faces%face_local(n,HI_SIDE,dir)
                             llo = faces%face_local(n,LO_SIDE,dir)
                             cell_val_flcl(lhi,1,nm) = min(cell_val_flcl(lhi,1,nm), &
                                  get_value(llo,nm))
                             cell_val_flcl(lhi,2,nm) = max(cell_val_flcl(lhi,2,nm), &
                                  get_value(llo,nm))
                          enddo ! n

                          !                   LO_SIDE
                          do n = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                             llo = faces%face_local(n,LO_SIDE,dir)
                             lhi = faces%face_local(n,HI_SIDE,dir)
                             cell_val_flcl(llo,1,nm) = min(cell_val_flcl(llo,1,nm), &
                                  get_value(lhi,nm))
                             cell_val_flcl(llo,2,nm) = max(cell_val_flcl(llo,2,nm), &
                                  get_value(lhi,nm))
                          enddo ! n

                       else if (faces%face_id(loop,dir).eq.2) then

                          if (kode(dir,nm).eq.-1 .or.                           &
                               & (kode(dir,nm).eq.-10 .and. faces%face_flag(loop,dir).eq.2) .or. &
                               & (kode(dir,nm).eq.-20 .and. faces%face_flag(loop,dir).eq.0)) then
                             do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                                llo = faces%face_local(n,LO_SIDE,dir)
                                cell_val_flcl(llo,1,nm) = min(cell_val_flcl(llo,1,nm), &
                                     -get_value(llo,nm))
                                cell_val_flcl(llo,2,nm) = max(cell_val_flcl(llo,2,nm), &
                                     -get_value(llo,nm))
                             enddo ! n
                          endif

                       else if (faces%face_id(loop,dir).eq.1) then

                          if (kode(dir,nm).eq.-1 .or. &
                               & (kode(dir,nm).eq.-10 .and. faces%face_flag(loop,dir).eq.2) .or. &
                               & (kode(dir,nm).eq.-20 .and. faces%face_flag(loop,dir).eq.0)) then
                             do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                                lhi = faces%face_local(n,HI_SIDE,dir)
                                cell_val_flcl(lhi,1,nm) = min(cell_val_flcl(lhi,1,nm), &
                                     -get_value(lhi,nm))
                                cell_val_flcl(lhi,2,nm) = max(cell_val_flcl(lhi,2,nm), &
                                     -get_value(lhi,nm))
                             enddo ! n
                          endif
                       endif ! face_id
                    enddo ! loop
                 enddo ! dir
              endif
           enddo ! nm
        endif ! (itrmax .ge. 1)

        TIMERSET(.false.,inside_com1)

      end associate

    end subroutine inside_com1_split
    ! ------------------------------------------------------------------------------
    subroutine inside_com2_split(mesh, itr, dir, numvec, deriv, deriv_mm, &
         & cell_value_hilo, invalue, value_cloned)

      use mesh_types, only : mesh_t

      ! -----  caling arguments
      type(mesh_t), intent(in)  :: mesh
      integer,      intent(in)  :: itr, dir
      integer,      intent(in)  :: numvec

      real(REAL64), intent(in)    :: deriv(:,:,:)
      real(REAL64), intent(inout) :: deriv_mm(:,:)
      real(REAL64), intent(inout) :: cell_value_hilo(:,:,:)
      real(REAL64), intent(in), dimension(:,:)           :: invalue
      real(REAL64), intent(in), optional, dimension(:,:) :: value_cloned

      !        local variables
      integer      :: nt, l, nm

      associate (cells => mesh%cells)

        TIMERSET(.true., inside_com2)

        ! ...  These lines will change when we finish threading the routine:

        if (itr .eq. 1) then
           do nm   = 1,numvec
              do nt = 1,mesh%levels%allnumtop
                 l   = mesh%levels%alltop(nt)
                 deriv_mm(l,nm) = deriv(l,dir,nm)
              enddo ! nt
           enddo ! nm
        endif ! itr

        if (dir .eq. 1) then
           do nm   = 1,numvec
              do l  = 1, mesh%cells%numcell_clone
                 cell_value_hilo(l,1,nm) = get_value(l,nm)
                 cell_value_hilo(l,2,nm) = get_value(l,nm)
              enddo ! nt
           enddo ! nm
        endif

        do nm    = 1,numvec
           do nt = 1,mesh%levels%allnumtop
              l   = mesh%levels%alltop(nt)
              cell_value_hilo(l,1,nm) = get_value(l,nm)- &
                   deriv(l,dir,nm)*cells%cell_half_lo(l,dir)
              cell_value_hilo(l,2,nm) = get_value(l,nm)+&
                   deriv(l,dir,nm)*cells%cell_half_hi(l,dir)
           enddo ! nt
        enddo ! nm

        TIMERSET(.false., inside_com2)

      end associate

    end subroutine inside_com2_split
    ! ------------------------------------------------------------------------------
    subroutine inside_com3_split(sim, mesh, itr, dir, numitr, numvec, &
         & lv_weight, kode, deriv, deriv_mm, cell_value_hilo, invalue, &
         & cell_deriv_hilo_all_dir, do_special, do_pressure, frac_core, core, &
         & gradp, intopt, deriv_weight, value_cloned, faceval )

      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t
      use mesh_state_types,      only : mesh_state_frac_core_t, mesh_state_core_t
      use gradient_types,        only : gradient_prop_t
      use interface_types,       only : interface_option_t

      ! non-optional scalars and derived types
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      integer,      intent(in)     :: itr, dir, numitr, numvec
      real(REAL64), intent(in)     :: lv_weight

      ! non-optional arrays
      integer,      intent(in)     :: kode(:,:)
      real(REAL64), intent(out)    :: deriv(:,:,:)
      real(REAL64), intent(inout) :: deriv_mm(:,:)
      real(REAL64), intent(inout) :: cell_value_hilo(:,:,:)
      real(REAL64), intent(in), dimension(:,:) :: invalue
      real(REAL64), intent(inout) :: cell_deriv_hilo_all_dir(:,:,:,:) !JV

      ! optional scalars and derived types
      logical,      intent(in), optional       :: do_special
      logical,      intent(in), optional       :: do_pressure
      type(mesh_state_frac_core_t), optional, intent(in) :: frac_core
      type(mesh_state_core_t), optional, intent(in) :: core
      type(gradient_prop_t), optional, intent(in) :: gradp
      type(interface_option_t), optional, intent(in) :: intopt

      ! optional arrays
      real(REAL64), intent(in), optional :: deriv_weight(:)
      real(REAL64), intent(in), optional, dimension(:,:) :: value_cloned
      !       for inflow bndy cond'ns
      real(REAL64), intent(in), optional     :: faceval(:,:)

      ! ----- local variables

      integer      :: nm, nt, l, loop, n, llo, lhi, dnr
      integer      :: lm, lmo, lmi
      real(REAL64) :: face_value
      real(REAL64) :: cell_val_mnmx_hilo(mesh%cells%numcell_clone,4)
      real(REAL64) :: cell_deriv_hilo(mesh%cells%numcell_clone,2)
      real(REAL64) :: vof_weight(mesh%cells%numcell_clone)

      associate (cells => mesh%cells, &
           faces => mesh%faces)

        TIMERSET(.true., inside_com3)

        do nm = 1,numvec
           ! .....    calculate the area weighted average face values

           if (mesh%levels%allnumtop.gt.0) then

              ! A
              ! .....    Define all the min's and max values to have the ones that will not change
              do l  = 1,mesh%cells%numcell_clone
                 cell_val_mnmx_hilo(l,1) = get_value(l,nm)
                 cell_val_mnmx_hilo(l,3) = get_value(l,nm)
                 cell_val_mnmx_hilo(l,2) = get_value(l,nm)
                 cell_val_mnmx_hilo(l,4) = get_value(l,nm)
              enddo

              do loop = 1,faces%face_num(dir)

                 if (faces%face_id(loop,dir) .gt. 2) then
                    !                 Do the HI_SIDE first
                    do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                       lmi = faces%face_local(n,HI_SIDE,dir)
                       cell_val_mnmx_hilo(lmi,1) =  huge(ONE)
                       cell_val_mnmx_hilo(lmi,2) = -huge(ONE)
                    enddo ! n
                    !                 Do the LO_SIDE
                    do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                       lmo = faces%face_local(n,LO_SIDE,dir)
                       cell_val_mnmx_hilo(lmo,3) =  huge(ONE)
                       cell_val_mnmx_hilo(lmo,4) = -huge(ONE)
                    enddo ! n
                 else if (faces%face_id(loop,dir) .eq. 2) then
                    do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                       lmo = faces%face_local(n,LO_SIDE,dir)
                       cell_val_mnmx_hilo(lmo,3) =  huge(ONE)
                       cell_val_mnmx_hilo(lmo,4) = -huge(ONE)
                    enddo ! n
                 else if (faces%face_id(loop,dir) .eq. 1) then
                    do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                       lmi = faces%face_local(n,HI_SIDE,dir)
                       cell_val_mnmx_hilo(lmi,1) =  huge(ONE)
                       cell_val_mnmx_hilo(lmi,2) = -huge(ONE)
                    enddo ! n
                 endif ! face_id
              enddo ! loop

              ! B
              ! .....    calculate the area weighted average face values
              do loop = 1,faces%face_num(dir)

                 if (faces%face_id(loop,dir) .gt. 2) then

                    if (present(do_special)) then
                       ! Do the HI_SIDE part
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          lhi = faces%face_local(n,HI_SIDE,dir)
                          llo = faces%face_local(n,LO_SIDE,dir)
                          face_value = ZERO

                          if (core%rho(llo).gt.ZERO .or. core%rho(lhi).gt.ZERO) then
                             if (do_pressure .and. cell_value_hilo(llo,1,nm)       &
                                  & * cell_value_hilo(lhi,2,nm).le.ZERO) then
                                ! .....  this coding addresses the hot spot problem
                                face_value = (cells%cell_half_lo(lhi,dir)*core%rho(lhi) *  &
                                     & cell_value_hilo(llo,2,nm)         &
                                     & +  cells%cell_half_hi(llo,dir)*core%rho(llo) *  &
                                     & cell_value_hilo(lhi,1,nm))        &
                                     & / (cells%cell_half_lo(lhi,dir)*core%rho(lhi)    &
                                     & +  cells%cell_half_hi(llo,dir)*core%rho(llo))
                             else
                                face_value = (cells%cell_half_lo(lhi,dir) *           &
                                     & cell_value_hilo(llo,2,nm)         &
                                     & +  cells%cell_half_hi(llo,dir) *           &
                                     & cell_value_hilo(lhi,1,nm))        &
                                     & / (cells%cell_half_lo(lhi,dir)             &
                                     & +  cells%cell_half_hi(llo,dir))
                             endif ! do_pressure ...
                          endif ! core%rho
                          cell_val_mnmx_hilo(lhi,1) =                          &
                               & min(cell_val_mnmx_hilo(lhi,1), face_value)
                          cell_val_mnmx_hilo(lhi,2) =                          &
                               & max(cell_val_mnmx_hilo(lhi,2), face_value)
                       enddo ! n

                       ! do the LO_SIDE part
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          llo = faces%face_local(n,LO_SIDE,dir)
                          lhi = faces%face_local(n,HI_SIDE,dir)
                          face_value = ZERO

                          if (core%rho(llo).gt.ZERO .or. core%rho(lhi).gt.ZERO) then
                             if (do_pressure .and. cell_value_hilo(llo,1,nm)       &
                                  & *cell_value_hilo(lhi,2,nm).le.ZERO) then
                                ! .....  this coding addresses the hot spot problem
                                face_value = (cells%cell_half_lo(lhi,dir)*core%rho(lhi) *  &
                                     & cell_value_hilo(llo,2,nm)         &
                                     & +  cells%cell_half_hi(llo,dir)*core%rho(llo) *  &
                                     & cell_value_hilo(lhi,1,nm))        &
                                     & / (cells%cell_half_lo(lhi,dir)*core%rho(lhi)    &
                                     & +  cells%cell_half_hi(llo,dir)*core%rho(llo))
                             else
                                face_value = (cells%cell_half_lo(lhi,dir) *           &
                                     & cell_value_hilo(llo,2,nm)         &
                                     & +  cells%cell_half_hi(llo,dir) *           &
                                     & cell_value_hilo(lhi,1,nm))        &
                                     & / (cells%cell_half_lo(lhi,dir)             &
                                     & +  cells%cell_half_hi(llo,dir))
                             endif
                          endif
                          cell_val_mnmx_hilo(llo,3) =                          &
                               & min(cell_val_mnmx_hilo(llo,3), face_value)
                          cell_val_mnmx_hilo(llo,4) =                          &
                               & max(cell_val_mnmx_hilo(llo,4), face_value)
                       enddo ! n

                    else !not present(do_special)
                       !       Do the HI_SIDE first
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          lhi = faces%face_local(n,HI_SIDE,dir)
                          llo = faces%face_local(n,LO_SIDE,dir)
                          face_value = (cells%cell_half_lo(lhi,dir)*cell_value_hilo(llo,2,nm)  &
                               & +cells%cell_half_hi(llo,dir)*cell_value_hilo(lhi,1,nm)) &
                               & /(cells%cell_half_lo(lhi,dir)+cells%cell_half_hi(llo,dir))
                          cell_val_mnmx_hilo(lhi,1) = min(cell_val_mnmx_hilo(lhi,1), face_value)
                          cell_val_mnmx_hilo(lhi,2) = max(cell_val_mnmx_hilo(lhi,2), face_value)
                       enddo ! n

                       !       Do the LO_SIDE
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          llo = faces%face_local(n,LO_SIDE,dir)
                          lhi = faces%face_local(n,HI_SIDE,dir)
                          face_value = (cells%cell_half_lo(lhi,dir)*cell_value_hilo(llo,2,nm)  &
                               & +cells%cell_half_hi(llo,dir)*cell_value_hilo(lhi,1,nm)) &
                               & /(cells%cell_half_lo(lhi,dir)+cells%cell_half_hi(llo,dir))
                          cell_val_mnmx_hilo(llo,3) = min(cell_val_mnmx_hilo(llo,3), face_value)
                          cell_val_mnmx_hilo(llo,4) = max(cell_val_mnmx_hilo(llo,4), face_value)
                       enddo ! n
                    endif ! present(do_special)

                 else if (faces%face_id(loop,dir) .eq. 2) then

                    if(kode(dir,nm).eq.-2)then
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          lmo = faces%face_local(n,LO_SIDE,dir)
                          face_value = faceval(n,dir)
                          cell_val_mnmx_hilo(lmo,3) = min(cell_val_mnmx_hilo(lmo,3), face_value)
                          cell_val_mnmx_hilo(lmo,4) = max(cell_val_mnmx_hilo(lmo,4), face_value)
                       enddo ! n
                    else
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          llo = faces%face_local(n,LO_SIDE,dir)
                          face_value = cell_value_hilo(llo,2,nm)
                          cell_val_mnmx_hilo(llo,3) = min(cell_val_mnmx_hilo(llo,3), face_value)
                          cell_val_mnmx_hilo(llo,4) = max(cell_val_mnmx_hilo(llo,4), face_value)
                       enddo ! n
                    endif

                 else if (faces%face_id(loop,dir) .eq. 1) then

                    if(kode(dir,nm).eq.-2)then
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          lmi = faces%face_local(n,HI_SIDE,dir)
                          face_value = faceval(n,dir)
                          cell_val_mnmx_hilo(lmi,1) = min(cell_val_mnmx_hilo(lmi,1), face_value)
                          cell_val_mnmx_hilo(lmi,2) = max(cell_val_mnmx_hilo(lmi,2), face_value)
                       enddo ! n
                    else
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          lhi = faces%face_local(n,HI_SIDE,dir)
                          face_value = cell_value_hilo(lhi,1,nm)
                          cell_val_mnmx_hilo(lhi,1) = min(cell_val_mnmx_hilo(lhi,1), face_value)
                          cell_val_mnmx_hilo(lhi,2) = max(cell_val_mnmx_hilo(lhi,2), face_value)
                       enddo ! n
                    endif
                 endif ! face_id

              enddo ! loop

              ! C
              do nt = 1,mesh%levels%allnumtop
                 lm  = mesh%levels%alltop(nt)
                 cell_deriv_hilo(lm,1:2) = ZERO
              enddo ! l

              do nt = 1,mesh%levels%allnumtop
                 l   = mesh%levels%alltop(nt)
                 if(cell_val_mnmx_hilo(l,1).gt. &
                      get_value(l,nm) &
                      )then
                    cell_deriv_hilo(l,1) = (get_value(l,nm) &
                         -cell_val_mnmx_hilo(l,1))  &
                         & / cells%cell_half_lo(l,dir)
                 else if(cell_val_mnmx_hilo(l,2).lt. &
                      get_value(l,nm) &
                      )then
                    cell_deriv_hilo(l,1) = (get_value(l,nm) &
                         -cell_val_mnmx_hilo(l,2))  &
                         & / cells%cell_half_lo(l,dir)
                 else
                    cell_deriv_hilo(l,1) = ZERO
                 endif

                 if(cell_val_mnmx_hilo(l,3).gt. &
                      get_value(l,nm) &
                      )then
                    cell_deriv_hilo(l,2) = (cell_val_mnmx_hilo(l,3)-&
                         get_value(l,nm) )  &
                         & / cells%cell_half_hi(l,dir)
                 else if(cell_val_mnmx_hilo(l,4).lt.&
                      get_value(l,nm) )then
                    cell_deriv_hilo(l,2) = (cell_val_mnmx_hilo(l,4)- &
                         get_value(l,nm) &
                         )  &
                         & / cells%cell_half_hi(l,dir)
                 else
                    cell_deriv_hilo(l,2) = ZERO
                 endif
              enddo ! nt

              ! D
              ! .....       restrict the derivatives at boundaries
              select case (kode(dir,nm))
              case (-1, -10, -20)

                 ! .....     quantities like velocity have to be zero at a boundary

                 do loop = 1,faces%face_num(dir)

                    if (kode(dir,nm).eq.-1 .or. &
                         & (kode(dir,nm).eq.-10 .and. faces%face_flag(loop,dir).eq.2) .or. &
                         & (kode(dir,nm).eq.-20 .and. faces%face_flag(loop,dir).eq.0)) then

                       if (faces%face_id(loop,dir) .eq. 1) then
                          do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                             lhi = faces%face_local(n,HI_SIDE,dir)
                             cell_deriv_hilo(lhi,1) = get_value(lhi,nm)&
                                  /cells%cell_half_lo(lhi,dir)
                          enddo ! n

                       else if (faces%face_id(loop,dir) .eq. 2) then
                          do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                             llo = faces%face_local(n,LO_SIDE,dir)
                             cell_deriv_hilo(llo,2) = -get_value(llo,nm)&
                                  /cells%cell_half_hi(llo,dir)
                          enddo ! n
                       endif !face_id
                    endif ! kode
                 enddo ! loop

              case (1)

                 do loop = 1,faces%face_num(dir)
                    if (faces%face_id(loop,dir) .eq. 1) then
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          lhi = faces%face_local(n,HI_SIDE,dir)
                          cell_deriv_hilo(lhi,1) =                                 &
                               & min(cell_deriv_hilo(lhi,2), &
                               get_value(lhi,nm)&
                               /cells%cell_half_lo(lhi,dir))
                       enddo ! n
                    else if (faces%face_id(loop,dir) .eq. 2) then
                       do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                          llo = faces%face_local(n,LO_SIDE,dir)
                          cell_deriv_hilo(llo,2) =                                 &
                               & max(cell_deriv_hilo(llo,1), &
                               -get_value(llo,nm)&
                               /cells%cell_half_hi(llo,dir))
                       enddo ! n
                    endif !face_id
                 enddo ! loop

              case (2, 3)

                 ! .....         gradients of pressure are continuous at gravitational boundaries

                 do loop = 1,faces%face_num(dir)
                    if (kode(dir,nm).eq.2 .or. &
                         & (kode(dir,nm).eq.3 .and. faces%face_flag(loop,dir).eq.1)) then

                       if (faces%face_id(loop,dir) .eq. 1) then
                          do n  = faces%face_lo(loop,dir),faces%face_hi(loop,dir)
                             lmi = faces%face_local(n,HI_SIDE,dir)
                             cell_deriv_hilo(lmi,1) = cell_deriv_hilo(lmi,2)
                          enddo ! n
                       else if (faces%face_id(loop,dir) .eq. 2) then
                          do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                             lmo = faces%face_local(n,LO_SIDE,dir)
                             cell_deriv_hilo(lmo,2) = cell_deriv_hilo(lmo,1)
                          enddo ! n
                       endif !face_id
                    endif !kode
                 enddo ! loop

              end select

              ! E
              ! .....       
              if (limit_slope .and. itr.eq.1) then
                 do nt = 1,mesh%levels%allnumtop
                    l   = mesh%levels%alltop(nt)
                    deriv_mm(l,nm) = max(ZERO,min(cell_deriv_hilo(l,1),cell_deriv_hilo(l,2))) &
                         & +min(ZERO,max(cell_deriv_hilo(l,1),cell_deriv_hilo(l,2)))
                 enddo ! nt
              endif

              ! .....       calculate the cell derivative
              ! F
              select case (method)
              case (INT_MM)  
                 do nt = 1,mesh%levels%allnumtop
                    l   = mesh%levels%alltop(nt)
                    deriv(l,dir,nm) = max(ZERO,min(cell_deriv_hilo(l,1),cell_deriv_hilo(l,2))) &
                         & +min(ZERO,max(cell_deriv_hilo(l,1),cell_deriv_hilo(l,2)))
                 enddo ! nt

              case (INT_EMM)  

                 do nt = 1,mesh%levels%allnumtop
                    l   = mesh%levels%alltop(nt)
                    deriv(l,dir,nm) = median(cell_deriv_hilo(l,1), &
                         & cell_deriv_hilo(l,2), &
                         & -(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))
                 enddo ! nt

              case (INT_LV)

                 if( numitr.eq.7 .or. numitr.eq.8 ) then !JV store left and right derivatives
                    do nt = 1,mesh%levels%allnumtop
                       l   = mesh%levels%alltop(nt)
                       cell_deriv_hilo_all_dir(dir,l,1:2,nm) = cell_deriv_hilo(l,1:2)
                       deriv(l,dir,nm) = (max(ZERO, min(lv_weight*cell_deriv_hilo(l,1), &
                            & lv_weight*cell_deriv_hilo(l,2), &
                            & HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))) &
                            & +min(ZERO, max(lv_weight*cell_deriv_hilo(l,1), &
                            & lv_weight*cell_deriv_hilo(l,2), &
                            & HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))))
                    enddo ! nt
                 else 
                    do nt = 1,mesh%levels%allnumtop
                       l   = mesh%levels%alltop(nt)
                       deriv(l,dir,nm) = (max(ZERO, min(lv_weight*cell_deriv_hilo(l,1), &
                            & lv_weight*cell_deriv_hilo(l,2), &
                            & HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))) &
                            & +min(ZERO, max(lv_weight*cell_deriv_hilo(l,1), &
                            & lv_weight*cell_deriv_hilo(l,2), &
                            & HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2)))))
                    enddo ! nt
                 endif

              case (INT_NL) ! --- turn off limiters

                 do nt = 1,mesh%levels%allnumtop
                    l   = mesh%levels%alltop(nt)
                    deriv(l,dir,nm) = HALF*(cell_deriv_hilo(l,1)+cell_deriv_hilo(l,2))
                 enddo ! nt

              end select

           endif !mesh%levels%allnumtop
        enddo ! nm

        ! G
        ! 
        if (itr .eq. itrmax) then
           if (present(gradp)) then
              if (limit_slope .and. present(deriv_weight) .and. gradp%shock_detector.ge.1) then
                 do nm   = 1,numvec
                    do nt = 1,mesh%levels%allnumtop
                       l   = mesh%levels%alltop(nt)
                       deriv(l,dir,nm) = deriv_mm(l,nm) +          &
                            & deriv_weight(l)*(deriv(l,dir,nm)-deriv_mm(l,nm))
                    enddo ! nt
                 enddo ! nm
              endif
           endif
        endif

        ! H
        if (dir .ge. sim%numdim) then
           if (present(intopt)) then
              ! ..... limit the slope to ZERO at interface fronts
              if (intopt%interface_option.ge.3) then

                 ! ..... Get where are the interface cells
                 vof_weight(1:mesh%cells%numcell_clone) = ONE
                 if (numitr.ge.0) call modify_weight_tempo(mesh, frac_core, intopt, vof_weight)
                 !-mf mar07
                 do nm    = 1,numvec
                    do dnr = 1,sim%numdim
                       if ((kode(dnr,nm).eq.-1.or.kode(dnr,nm).eq.-10 .or.                  &
                            & kode(dnr,nm).eq.30.or.kode(dnr,nm).eq.2).and. &
                            & .not. intopt%flatten_interface_vel) then
                          !-vel, pres, grav do nothing
                       else
                          do n  = 1,mesh%levels%allnumtop
                             l  = mesh%levels%alltop(n)
                             deriv(l,dnr,nm) = deriv(l,dnr,nm)*vof_weight(l)
                          enddo ! nt
                       endif
                    enddo
                 enddo
              endif
           endif
        endif

        TIMERSET(.false., inside_com3)

      end associate

    end subroutine inside_com3_split
    ! ------------------------------------------------------------------------------
    subroutine inside_com3A_split(sim, mesh, dir, nm, cell_val_mnmx_hilo, invalue, &
         & value_cloned )

      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t

      ! non-optional scalars and derived types
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      integer,      intent(in)     :: dir
      integer,      intent(in)     :: nm

      ! non-optional arrays
      real(REAL64), intent(inout), dimension(:,:) :: cell_val_mnmx_hilo !(mesh%cells%numcell_clone,4)
      real(REAL64), intent(in), dimension(:,:) :: invalue

      ! optional arrays
      real(REAL64), intent(in), optional, dimension(:,:) :: value_cloned

      ! ----- local variables
      integer      :: l, loop, lmo, lmi, n

      associate (cells => mesh%cells, &
           & faces => mesh%faces)

        TIMERSET(.true., inside_com3a)

        ! .....    Define all the min's and max values to have the ones that will not change
        do l  = 1,mesh%cells%numcell_clone
           cell_val_mnmx_hilo(l,1) = get_value(l,nm)
           cell_val_mnmx_hilo(l,3) = get_value(l,nm)
           cell_val_mnmx_hilo(l,2) = get_value(l,nm)
           cell_val_mnmx_hilo(l,4) = get_value(l,nm)
        enddo

        do loop = 1,faces%face_num(dir)

           if (faces%face_id(loop,dir) .gt. 2) then
              !                 Do the HI_SIDE first
              do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                 lmi = faces%face_local(n,HI_SIDE,dir)
                 cell_val_mnmx_hilo(lmi,1) =  huge(ONE)
                 cell_val_mnmx_hilo(lmi,2) = -huge(ONE)
              enddo ! n
              !                 Do the LO_SIDE
              do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                 lmo = faces%face_local(n,LO_SIDE,dir)
                 cell_val_mnmx_hilo(lmo,3) =  huge(ONE)
                 cell_val_mnmx_hilo(lmo,4) = -huge(ONE)
              enddo ! n
           else if (faces%face_id(loop,dir) .eq. 2) then
              do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                 lmo = faces%face_local(n,LO_SIDE,dir)
                 cell_val_mnmx_hilo(lmo,3) =  huge(ONE)
                 cell_val_mnmx_hilo(lmo,4) = -huge(ONE)
              enddo ! n
           else if (faces%face_id(loop,dir) .eq. 1) then
              do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                 lmi = faces%face_local(n,HI_SIDE,dir)
                 cell_val_mnmx_hilo(lmi,1) =  huge(ONE)
                 cell_val_mnmx_hilo(lmi,2) = -huge(ONE)
              enddo ! n
           endif ! face_id
        enddo ! loop

        TIMERSET(.false., inside_com3a)

      end associate

    end subroutine inside_com3A_split
    ! ------------------------------------------------------------------------------
    subroutine inside_com3C_split(sim, mesh, dir, nm, cell_deriv_hilo, &
         & cell_val_mnmx_hilo, invalue, value_cloned )

      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t

      ! non-optional scalars and derived types
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      integer,      intent(in)     :: dir, nm

      ! non-optional arrays
      real(REAL64), intent(inout), dimension(:,:) :: cell_deriv_hilo !(mesh%cells%numcell_clone,2)
      real(REAL64), intent(inout), dimension(:,:) :: cell_val_mnmx_hilo !(mesh%cells%numcell_clone,4)
      real(REAL64), intent(in), dimension(:,:) :: invalue

      ! optional arrays
      real(REAL64), intent(in), optional, dimension(:,:) :: value_cloned

      ! ----- local variables
      integer :: nt, l, lm

      associate (cells => mesh%cells, &
           faces => mesh%faces)

        TIMERSET(.true., inside_com3c)

        do nt = 1,mesh%levels%allnumtop
           lm  = mesh%levels%alltop(nt)
           cell_deriv_hilo(lm,1) = ZERO
           cell_deriv_hilo(lm,2) = ZERO
        enddo ! l

        do nt = 1,mesh%levels%allnumtop
           l   = mesh%levels%alltop(nt)
           if(cell_val_mnmx_hilo(l,1).gt.get_value(l,nm))then
              cell_deriv_hilo(l,1) = (get_value(l,nm) &
                   -cell_val_mnmx_hilo(l,1))  &
                   & / cells%cell_half_lo(l,dir)
           else if(cell_val_mnmx_hilo(l,2).lt.get_value(l,nm))then
              cell_deriv_hilo(l,1) = (get_value(l,nm) &
                   -cell_val_mnmx_hilo(l,2))  &
                   & / cells%cell_half_lo(l,dir)
           else
              cell_deriv_hilo(l,1) = ZERO
           endif

           if(cell_val_mnmx_hilo(l,3).gt.get_value(l,nm))then
              cell_deriv_hilo(l,2) = (cell_val_mnmx_hilo(l,3)- &
                   get_value(l,nm) &
                   )  &
                   & / cells%cell_half_hi(l,dir)
           else if(cell_val_mnmx_hilo(l,4).lt.get_value(l,nm))then
              cell_deriv_hilo(l,2) = (cell_val_mnmx_hilo(l,4)- &
                   get_value(l,nm) &
                   )  &
                   & / cells%cell_half_hi(l,dir)
           else
              cell_deriv_hilo(l,2) = ZERO
           endif
        enddo ! nt

        TIMERSET(.false., inside_com3c)

      end associate

    end subroutine inside_com3C_split
    ! ------------------------------------------------------------------------------
    subroutine inside_com3D_split(sim, mesh, dir, nm, kode, cell_deriv_hilo, invalue, &
         & value_cloned )

      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t

      ! non-optional scalars and derived types
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      integer,      intent(in)     :: dir, nm

      ! non-optional arrays
      integer,      intent(in)     :: kode(:,:)
      real(REAL64), intent(inout), dimension(:,:) :: cell_deriv_hilo !(mesh%cells%numcell_clone,2)
      real(REAL64), intent(in), dimension(:,:) :: invalue

      ! optional arrays
      real(REAL64), intent(in), optional, dimension(:,:) :: value_cloned

      ! ----- local variables
      integer      :: loop, n, llo, lhi, lmo, lmi

      associate (cells => mesh%cells, &
           faces => mesh%faces)

        TIMERSET(.true., inside_com3d)

        ! .....       restrict the derivatives at boundaries
        select case (kode(dir,nm))
        case (-1, -10, -20)

           ! .....     quantities like velocity have to be zero at a boundary

           do loop = 1,faces%face_num(dir)

              if (kode(dir,nm).eq.-1 .or. &
                   & (kode(dir,nm).eq.-10 .and. faces%face_flag(loop,dir).eq.2) .or. &
                   & (kode(dir,nm).eq.-20 .and. faces%face_flag(loop,dir).eq.0)) then

                 if (faces%face_id(loop,dir) .eq. 1) then
                    do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                       lhi = faces%face_local(n,HI_SIDE,dir)
                       cell_deriv_hilo(lhi,1) = &
                            get_value(lhi,nm)/cells%cell_half_lo(lhi,dir)
                    enddo ! n

                 else if (faces%face_id(loop,dir) .eq. 2) then
                    do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                       llo = faces%face_local(n,LO_SIDE,dir)
                       cell_deriv_hilo(llo,2) = &
                            -get_value(llo,nm)/cells%cell_half_hi(llo,dir)
                    enddo ! n
                 endif !face_id
              endif ! kode
           enddo ! loop

        case (1)

           do loop = 1,faces%face_num(dir)
              if (faces%face_id(loop,dir) .eq. 1) then
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    lhi = faces%face_local(n,HI_SIDE,dir)
                    cell_deriv_hilo(lhi,1) =                                 &
                         & min(cell_deriv_hilo(lhi,2), &
                         get_value(lhi,nm)/&
                         cells%cell_half_lo(lhi,dir))
                 enddo ! n
              else if (faces%face_id(loop,dir) .eq. 2) then
                 do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                    llo = faces%face_local(n,LO_SIDE,dir)
                    cell_deriv_hilo(llo,2) =                                 &
                         & max(cell_deriv_hilo(llo,1), &
                         -get_value(llo,nm)/&
                         cells%cell_half_hi(llo,dir))
                 enddo ! n
              endif !face_id
           enddo ! loop

        case (2, 3)

           ! .....         gradients of pressure are continuous at gravitational boundaries

           do loop = 1,faces%face_num(dir)
              if (kode(dir,nm).eq.2 .or. &
                   & (kode(dir,nm).eq.3 .and. faces%face_flag(loop,dir).eq.1)) then

                 if (faces%face_id(loop,dir) .eq. 1) then
                    do n  = faces%face_lo(loop,dir),faces%face_hi(loop,dir)
                       lmi = faces%face_local(n,HI_SIDE,dir)
                       cell_deriv_hilo(lmi,1) = cell_deriv_hilo(lmi,2)
                    enddo ! n
                 else if (faces%face_id(loop,dir) .eq. 2) then
                    do n  = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
                       lmo = faces%face_local(n,LO_SIDE,dir)
                       cell_deriv_hilo(lmo,2) = cell_deriv_hilo(lmo,1)
                    enddo ! n
                 endif !face_id
              endif !kode
           enddo ! loop

        end select

        TIMERSET(.false., inside_com3d)

      end associate

    end subroutine inside_com3D_split
    ! ------------------------------------------------------------------------------
    subroutine inside_com4_split(sim, mesh, numvec, numitr, lv_weight, deriv, &
         & cell_val_flcl, cell_deriv_hilo_all_dir, invalue, value_cloned )

      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t

      ! -----  caling arguments
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in)      :: mesh
      integer,      intent(in)      :: numvec
      integer,      intent(in)      :: numitr                           !JV
      real(REAL64), intent(inout)   :: lv_weight

      real(REAL64), intent(inout) :: deriv(:,:,:)
      real(REAL64), intent(in), dimension(:,:,:) :: cell_val_flcl
      real(REAL64), intent(in) :: cell_deriv_hilo_all_dir(:,:,:,:) !JV
      real(REAL64), intent(in), dimension(:,:) :: invalue
      real(REAL64), intent(in), optional, dimension(:,:) :: value_cloned

      !cdsv todo value_c, array_id, material_id

      !    --- local variables
      real(REAL64) :: delta, scale
      integer      :: nt, l, dir, nm
      real(REAL64) :: normp !JV norm of the original gradient
      real(REAL64) :: lv_weight_loc !JV local lv_weight based on
      ! the direction of the original gradient (ranging from 1 to sqrt(3))
      associate ( cells => mesh%cells, &
           levs => mesh%levels)

        TIMERSET(.true., inside_com4)

        if (method.ne.INT_NL) then
           if (numitr .eq. 7 .or. numitr .eq. 8) then 
              do nm   = 1,numvec
                 do nt = 1,levs%allnumtop
                    l   = levs%alltop(nt)
                    normp = 0.0d0
                    do dir  = 1,sim%numdim
                       normp = normp + deriv(l,dir,nm)**2
                    enddo
                    normp = sqrt(normp)
                    if ( normp > 1.0d-13 ) then !JV
                       lv_weight_loc = 0.0d0
                       do dir  = 1,sim%numdim
                          lv_weight_loc = lv_weight_loc + ( abs( deriv(l,dir,nm) ) / normp )**3
                       enddo
                       lv_weight_loc = lv_weight_loc * lv_weight
                    else
                       lv_weight_loc = lv_weight
                    endif
                    if(numitr .eq. 8) lv_weight_loc = lv_weight !JV temporary test selection
                    do dir  = 1,sim%numdim
                       !JV to make sure that we do not increase the slopes
                       deriv(l,dir,nm) = sign( min( abs( deriv(l,dir,nm) ) , &
                            & abs(max(ZERO, min(lv_weight_loc*cell_deriv_hilo_all_dir(dir,l,1,nm), &
                            & lv_weight_loc*cell_deriv_hilo_all_dir(dir,l,2,nm), &
                            & HALF*(cell_deriv_hilo_all_dir(dir,l,1,nm)+cell_deriv_hilo_all_dir(dir,l,2,nm)))) &
                            & +min(ZERO, max(lv_weight_loc*cell_deriv_hilo_all_dir(dir,l,1,nm), &
                            & lv_weight_loc*cell_deriv_hilo_all_dir(dir,l,2,nm), &
                            & HALF*(cell_deriv_hilo_all_dir(dir,l,1,nm)+cell_deriv_hilo_all_dir(dir,l,2,nm))))) ) , &
                            & deriv(l,dir,nm) )
                    enddo !dir
                 enddo ! nt
              enddo ! nm
           else !JV original cell bound threatment
              do nm   = 1,numvec
                 do nt = 1,levs%allnumtop
                    l  = levs%alltop(nt)

                    delta   = ZERO
                    do dir  = 1,sim%numdim
                       delta = delta+abs(deriv(l,dir,nm))*cells%cell_half_lo(l,dir)
                    enddo ! dir
                    if(get_value(l,nm) &
                         +delta.gt.cell_val_flcl(l,2,nm) .or.   &
                         & get_value(l,nm) &
                         -delta.lt.cell_val_flcl(l,1,nm))then
                       scale = min(get_value(l,nm) &
                            -cell_val_flcl(l,1,nm),     &
                            & cell_val_flcl(l,2,nm)-&
                            get_value(l,nm))     &
                            & / delta
                       deriv(l,:sim%numdim,nm) = scale*deriv(l,:sim%numdim,nm)
                    endif
                 enddo ! nt
              enddo ! nm
           endif ! (numitr .eq. 7 .or. numitr .eq. 8)
        endif ! method

        TIMERSET(.false., inside_com4)

      end associate

    end subroutine inside_com4_split
    ! ------------------------------------------------------------------------------
    subroutine derivatives_common_split(sim, mesh, &
         frac_core, core, &
         gradp, intopt, &
         cell_dim, numitr, nvec, kode, &
         noslope_cell, deriv, do_fincom, faceval, &
         deriv_weight, invalue, value_cloned, do_special)
      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t
      use mesh_state_types, only : mesh_state_frac_core_t, mesh_state_core_t
      use gradient_types,         only : gradient_prop_t
      use interface_types,        only : interface_option_t
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      type(mesh_state_frac_core_t), intent(in) :: frac_core
      type(mesh_state_core_t), intent(in) :: core
      type(gradient_prop_t), intent(in) :: gradp
      type(interface_option_t), intent(in) :: intopt
      integer,     intent(in) :: cell_dim
      integer,     intent(in) :: numitr
      integer,     intent(in) :: nvec
      integer,     intent(in) :: kode(:,:)
      logical,     intent(in), allocatable :: noslope_cell(:)
      logical,     intent(in) :: do_fincom
      logical,     intent(in), optional :: do_special

      real(REAL64), intent(out),  dimension(:,:,:) :: deriv

      ! for inflow bndy cond'ns
      real(REAL64), intent(in), optional            :: faceval(:,:)
      real(REAL64), intent(in), optional            :: deriv_weight(:)

      real(REAL64), intent(in),    dimension(:,:)   :: invalue
      real(REAL64), intent(in), optional, dimension(:,:)   :: value_cloned
      
      call derivatives_common_internal_split(sim, mesh, cell_dim, numitr, nvec, kode, &
           noslope_cell, deriv, do_fincom, &
           invalue, &
           faceval=faceval, deriv_weight=deriv_weight, &
           value_cloned=value_cloned, do_special=do_special, &
           frac_core=frac_core, core=core, &
           gradp=gradp, intopt=intopt)
    end subroutine derivatives_common_split

    subroutine derivatives_common_simple_split(sim, mesh, &
         cell_dim, numitr, nvec, kode, &
         noslope_cell, deriv, do_fincom, faceval, &
         deriv_weight, invalue, value_cloned, do_special)
      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      integer,     intent(in) :: cell_dim
      integer,     intent(in) :: numitr
      integer,     intent(in) :: nvec
      integer,     intent(in) :: kode(:,:)
      logical,     intent(in), allocatable :: noslope_cell(:)
      logical,     intent(in) :: do_fincom
      logical,     intent(in), optional :: do_special

      real(REAL64), intent(out),  dimension(:,:,:) :: deriv

      ! for inflow bndy cond'ns
      real(REAL64), intent(in), optional            :: faceval(:,:)
      real(REAL64), intent(in), optional            :: deriv_weight(:)

      real(REAL64), intent(in),    dimension(:,:)   :: invalue
      real(REAL64), intent(in), optional, dimension(:,:)   :: value_cloned

      call derivatives_common_internal_split(sim, mesh, cell_dim, numitr, nvec, kode, &
           noslope_cell, deriv, do_fincom, &
           invalue, &
           faceval=faceval, deriv_weight=deriv_weight, &
           value_cloned=value_cloned, do_special=do_special)
    end subroutine derivatives_common_simple_split

    subroutine derivatives_common_internal_split(sim, mesh, &
         cell_dim, numitr, nvec, kode, &
         noslope_cell, deriv, do_fincom, &
         invalue, &
         faceval, deriv_weight, &
         value_cloned, do_special, &
         frac_core, core, &
         gradp, intopt)
      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t
      use mesh_state_types, only : mesh_state_frac_core_t, mesh_state_core_t
      use gradient_types,         only : gradient_prop_t
      use interface_types,        only : interface_option_t
      use clone_lib_module,       only : clone_get
      use timer_module  ,         only : timerset
      use util,                   only : global_error
      use var_wrapper_class,      only : var_wrapper, vw_set
#ifdef EAP_KOKKOS_GRADIENTS
      use flcl_mod, only: nd_array_t, to_nd_array
      use iso_c_binding
      use gradients_interfaces
#endif
      ! ------------------------------------------------------------------------------
      ! ----- Purpose
      !       This routine calculates the derivative for variable value for a
      !       variable with dimension (numcell_clone,nvec). When using OpenMP,
      !       it calculates the derivative of value, i.e. dvalue/ddirection,
      !       distributed into threads.

      !       In certain cases, this routine makes an approximation: It considers
      !       zero the derivatives that are very close to zero (TINY_DERIVATIVE).
      !       It will only save time if we have a large number of derivatives close to
      !       zero - abs(value(a)-value(b))/abs(value(a)+value(b)) .le. TINY_DERIVATIVE -
      !       it will make these derivatives = 0 without trying to calculate the exact
      !       value.
      !       The subroutine deriv_details specify how we are calculating the derivatives.
      ! ------------------------------------------------------------------------------

      ! ----- calling arguments
      class(sim_info_t), intent(in)                :: sim
      type(mesh_t), intent(in)                     :: mesh
      integer,     intent(in)                      :: cell_dim
      integer,     intent(in)                      :: numitr
      integer,     intent(in)                      :: nvec
      integer,     intent(in)                      :: kode(:,:)
      logical,     intent(in), allocatable         :: noslope_cell(:)
      logical,     intent(in)                      :: do_fincom
      logical,     intent(in),    optional         :: do_special
      type(mesh_state_frac_core_t), optional, intent(in) :: frac_core
      type(mesh_state_core_t), optional, intent(in) :: core
      type(gradient_prop_t), optional, intent(in)                 :: gradp
      type(interface_option_t), optional, intent(in)           :: intopt

      real(REAL64), intent(out),  dimension(:,:,:) :: deriv

      ! for inflow bndy cond'ns
      real(REAL64), intent(in), optional            :: faceval(:,:)
      real(REAL64), intent(in), optional            :: deriv_weight(:)

      real(REAL64), intent(in),    dimension(:,:)   :: invalue
      real(REAL64), intent(in), optional, dimension(:,:)   :: value_cloned


      ! ----- local variables

      logical      :: do_pressure

      integer      :: dir, itr, ivec, nm
      integer      :: lvofd
      integer      :: priv_suntop, suntop
      real(REAL64) :: lv_weight = ZERO
      type(var_wrapper) :: vw(1)

      real(REAL64) :: cell_val_flcl(mesh%cells%numcell_clone,2,nvec)
      real(REAL64) :: cell_value_hilo(mesh%cells%numcell_clone,2,nvec)
      real(REAL64) :: deriv_mm(mesh%cells%numcell_clone,nvec)
      real(REAL64) :: cell_val_mnmx_hilo(mesh%cells%numcell_clone,4)
      real(REAL64) :: cell_deriv_hilo(mesh%cells%numcell_clone,2)

      real(REAL64), allocatable :: cell_deriv_hilo_all_dir(:,:,:,:)    !JV left and right extimates of derivatives

#ifdef EAP_KOKKOS_GRADIENTS
      ! ----- interfacing variables

      logical(c_bool) :: ic_present_do_special
      logical(c_bool) :: ic_do_special
      logical(c_bool) :: ic_present_do_pressure
      logical(c_bool) :: ic_do_pressure
      logical(c_bool) :: ic_present_gradp
      logical(c_bool) :: ic_present_deriv_weight
      logical(c_bool) :: mype_is_iope
      logical(c_bool) :: ic_present_faceval
      logical(c_bool) :: ic_limit_slope
      type(nd_array_t) :: ic_faceval_wrapper
      type(nd_array_t) :: ic_deriv_weight_wrapper
      type(gradient_prop_t) :: ic_gradp

      ! ----- dummy variables
      real(REAL64) :: fake_faceval(2,2)
      real(REAL64) :: fake_deriv_weight(2)

#endif

      ! ----- code

      associate (cells => mesh%cells, &
           levs => mesh%levels, &
           faces => mesh%faces)

        TIMERSET_UNCOND(.true.,main)

        !JV allocate an array for left and right extimates of derivatives
        if( numitr.eq.7 .or. numitr.eq.8 ) allocate( cell_deriv_hilo_all_dir(sim%numdim,mesh%cells%numcell_clone,2,nvec) )

        call deriv_details(numitr, lv_weight, cell_dim, do_pressure,  &
             do_special=do_special)

#ifdef EAP_KOKKOS_GRADIENTS
        ! set up interfacing variables after deriv_details
        if( present(faceval) ) then
           ic_faceval_wrapper = to_nd_array(faceval)
           ic_present_faceval = .true.
        else
           ic_faceval_wrapper = to_nd_array(fake_faceval)
           ic_present_faceval = .false.
        end if

        ic_present_do_pressure = .true.
        if ( do_pressure ) then
           ic_do_pressure = .true.
        else
           ic_do_pressure = .false.
        end if

        if(present(do_special)) then
           ic_present_do_special = .true.
           ic_do_special = logical(do_special, c_bool)
        else
           ic_present_do_special = .false.
           ic_do_special = .false.
        end if

        if(present(gradp)) then
           ic_present_gradp = .true.
           ic_gradp = gradp
        else
           ic_present_gradp = .false.
           ! ic3_gradp will be uninitialized, don't use it
        end if

        if(present(deriv_weight)) then
           ic_present_deriv_weight = .true.
           ic_deriv_weight_wrapper = to_nd_array(deriv_weight)
        else
           ic_present_deriv_weight = .false.
           ic_deriv_weight_wrapper = to_nd_array(fake_deriv_weight)
        end if

        mype_is_iope = .false.
        if (mype.eq.iope) mype_is_iope = .true.

        ic_limit_slope = logical(limit_slope, c_bool)
#endif
        call inside_com1_split(sim, mesh, core, nvec, &
             & kode, deriv, cell_val_flcl, invalue=invalue, value_cloned=value_cloned)

        if (itrmax.gt.0) then
           !         calculate all, from now on:

#ifdef _DEBUG
           do ivec   = 1,nvec
              do dir = 1,sim%numdim
                 if(kode(dir,ivec).eq.-2)then
                    if(.not.present(faceval))then
                       call global_error(                                    &
                            'DERIVATIVES: optional argument missing for kode = -2')
                    endif

                 endif
              enddo
           enddo
#endif

           lvofd = 1
           if (intopt%interface_option .ge. 3) lvofd = mesh%cells%numcell_clone

           do dir    = 1,sim%numdim
              do itr  = 1,itrmax
                 TIMERSET_UNCOND(.true.,inner)
                 ! ......      Estimate values at edges of cells
#ifdef EAP_KOKKOS_GRADIENTS
                 call inside_com2_split_arrays(  itr, dir, nvec, &
                      & mesh%levels%allnumtop, &
                      & mesh%cells%numcell_clone, &
                      & mesh%cells%numcell, mype_is_iope, &
                      & to_nd_array(mesh%levels%alltop), &
                      & to_nd_array(deriv), &
                      & to_nd_array(deriv_mm), &
                      & to_nd_array(cell_value_hilo), &
                      & to_nd_array(invalue), &
                      & to_nd_array(value_cloned), &
                      & to_nd_array(cells%cell_half_lo), &
                      & to_nd_array(cells%cell_half_hi) )
#else
                 ! else ifdef EAP_KOKKOS_GRADIENTS
                 ! no EAP_KOKKOS_GRADIENTS - so we're running F only code
                 call inside_com2_split(mesh, itr, dir, nvec, deriv, deriv_mm, &
                      & cell_value_hilo, invalue=invalue, value_cloned=value_cloned)
#endif
                 !endif EAP_KOKKOS_GRADIENTS

                 !             Not necessary to communicate cell_value_hilo here because they have
                 !             been estimated in the ghost cells.
#ifdef EAP_KOKKOS_GRADIENTS
                 ! don't zero out cell_val_mnmx_hilo (initialized in inside_com3a)
                 ! don't zero out cell_deriv_hilo (zeroed out in inside_com3c)
                 do nm = 1,nvec
                    ! .....    calculate the area weighted average face values
                    if (mesh%levels%allnumtop.gt.0) then
                       call inside_com3a_split_arrays( mype_is_iope, dir, nm, cells%numcell,&
                            & cells%numcell_clone, huge(ONE), to_nd_array(faces%face_num), &
                            & to_nd_array(faces%face_id), to_nd_array(faces%face_lo), &
                            & to_nd_array(faces%face_hi), to_nd_array(faces%face_local), &
                            & to_nd_array(cell_val_mnmx_hilo), to_nd_array(invalue), &
                            & to_nd_array(value_cloned) )
                       call inside_com3b_arrays( mype_is_iope, ic_present_faceval, &
                            & ic_present_do_special, ic_do_special, ic_present_do_pressure, &
                            & ic_do_pressure, dir, nm, to_nd_array(faces%face_num), &
                            & to_nd_array(faces%face_id), to_nd_array(faces%face_lo), &
                            & to_nd_array(faces%face_hi), &
                            & to_nd_array(faces%face_local), to_nd_array(kode), &
                            & to_nd_array(cell_val_mnmx_hilo), to_nd_array(cell_value_hilo), &
                            & ic_faceval_wrapper, to_nd_array(core%rho), &
                            & to_nd_array(cells%cell_half_lo), &
                            & to_nd_array(cells%cell_half_hi) )
                       call inside_com3c_split_arrays( mype_is_iope, dir, nm, &
                            & levs%allnumtop, cells%numcell, to_nd_array(levs%alltop), &
                            & to_nd_array(cells%cell_half_lo), to_nd_array(cells%cell_half_hi), &
                            & to_nd_array(cell_val_mnmx_hilo), to_nd_array(cell_deriv_hilo), &
                            & to_nd_array(invalue), to_nd_array(value_cloned) )
                       call inside_com3d_split_arrays( mype_is_iope, dir, nm, cells%numcell, &
                            & to_nd_array(faces%face_num), to_nd_array(faces%face_id), &
                            & to_nd_array(faces%face_lo), to_nd_array(faces%face_hi), &
                            & to_nd_array(faces%face_local), to_nd_array(faces%face_flag), &
                            & to_nd_array(kode), to_nd_array(cells%cell_half_lo), &
                            & to_nd_array(cells%cell_half_hi), to_nd_array(cell_deriv_hilo), &
                            & to_nd_array(invalue), to_nd_array(value_cloned) )
                       call inside_com3e_arrays( mype_is_iope, ic_limit_slope, itr, nm, &
                            & levs%allnumtop, to_nd_array(levs%alltop), &
                            & to_nd_array(cell_deriv_hilo), to_nd_array(deriv_mm) )
                       call inside_com3f_arrays( mype_is_iope, dir, nm, numitr, &
                            & method, levs%allnumtop, lv_weight, &
                            & to_nd_array(levs%alltop), to_nd_array(deriv), &
                            & to_nd_array(cell_deriv_hilo), to_nd_array(cell_deriv_hilo_all_dir) )
                    end if
                 end do
                 call inside_com3g_arrays( mype_is_iope, ic_limit_slope, &
                      & ic_present_gradp, ic_present_deriv_weight, dir, nvec, &
                      & itr, itrmax, levs%allnumtop, ic_gradp, &
                      & to_nd_array(levs%alltop), to_nd_array(deriv), &
                      & to_nd_array(deriv_mm), ic_deriv_weight_wrapper )
                 call inside_com3h(sim, mesh, dir, numitr, nvec, &
                      & kode, deriv, frac_core=frac_core, intopt=intopt )
#else
                 ! else ifdef EAP_KOKKOS_GRADIENTS
                 ! no EAP_KOKKOS_GRADIENTS - so we're running only F code
                 ! don't zero out cell_val_mnmx_hilo (initialized in inside_com3a)
                 ! don't zero out cell_deriv_hilo (zeroed out in inside_com3c)
                 do nm = 1,nvec
                    ! .....    calculate the area weighted average face values
                    if (mesh%levels%allnumtop.gt.0) then
                       call inside_com3A_split(sim, mesh, dir, nm, cell_val_mnmx_hilo, &
                            & invalue, value_cloned=value_cloned )
                       call inside_com3b(sim, mesh, dir, nm, cell_val_mnmx_hilo, &
                            & kode, cell_value_hilo, do_special=do_special, &
                            & do_pressure=do_pressure, core=core, faceval=faceval )
                       call inside_com3C_split(sim, mesh, dir, nm, cell_deriv_hilo, &
                            & cell_val_mnmx_hilo, invalue, value_cloned=value_cloned )
                       call inside_com3D_split(sim, mesh, dir, nm, kode, cell_deriv_hilo, &
                            & invalue, value_cloned=value_cloned )
                       call inside_com3e(sim, mesh, itr, nm, limit_slope, &
                            & cell_deriv_hilo, deriv_mm )
                       call inside_com3f(sim, mesh, dir, nm, numitr, method, &
                            & lv_weight, deriv, cell_deriv_hilo, &
                            & cell_deriv_hilo_all_dir )
                    end if
                 end do
                 call inside_com3g(mesh, dir, nvec, itr, itrmax, limit_slope, &
                      & deriv, deriv_mm, gradp=gradp, deriv_weight=deriv_weight )
                 call inside_com3h(sim, mesh, dir, numitr, nvec, &
                      & kode, deriv, frac_core=frac_core, intopt=intopt )
#endif
                 ! endif EAP_KOKKOS_GRADIENTS
                 TIMERSET_UNCOND(.false.,inner)
                 if (itr.lt.itrmax) then
                    TIMERSET(.true., clone_get_1)
                    call clone_get(deriv(:cells%numcell_clone,dir,:nvec),cells%numcell_clone,nvec)
                    TIMERSET(.false., clone_get_1)
                 endif

              enddo ! itr
           enddo ! dir

           call inside_com4_split(sim, mesh, nvec, numitr, lv_weight, deriv, &
                & cell_val_flcl, cell_deriv_hilo_all_dir, invalue, &
                & value_cloned=value_cloned )

           ! ....    deallocate variables allocated before com1:

        endif ! itrmax

        !       Do the communications necessary before calling hydro_face_nocomm
        !       Communicate both value and deriv before calling hydro_face_nocomm

        if (do_fincom) then
           !           We do not need to communicate "value" here. It has been already
           !           communicated in the beginning of the routine.
           TIMERSET(.true., clone_get_2)
           call vw_set(vw(1), deriv, cells%numcell, sim%numdim, nvec)
           call clone_get(vw, 1)
           TIMERSET(.false., clone_get_2)
        endif

        if( numitr.eq.7 .or. numitr.eq.8 ) deallocate( cell_deriv_hilo_all_dir ) !JV

        !SS TIMERSET(.true., post_modify_slope)
        !SS if (allocated(noslope_cell)) &
        !SS      & call post_modify_slope(cell_dim,sim%numdim,nvec,noslope_cell,deriv)
        !SS TIMERSET(.false., post_modify_slope)

        TIMERSET_UNCOND(.false.,main)

      end associate

    contains
      ! ------------------------------------------------------------------------------
      subroutine deriv_details(numitr, lv_weight, cell_dim, &
           & do_pressure, do_special)

        use util       , only : global_error

        integer,      intent(in)     :: cell_dim
        integer,      intent(in)     :: numitr
        logical,      intent(out)    :: do_pressure
        real(REAL64), intent(inout)  :: lv_weight

        logical,      intent(in), optional :: do_special

        logical     :: dimerror = .false.

        !       numitr = 0 : no derivatives (1st order accurate)
        !              = 1 : MM
        !              = 2 : iterated MM
        !              = 3 : extended MM
        !              = 4 : LV
        !              = 5 : NL
        !              = 6 : modified LV

        do_pressure = .false.
        if (present(do_special)) then
           do_pressure = do_special
        endif

        nitr = numitr
        if (nitr.lt.0) nitr = -nitr

        ! ..... select method
        method = NO_DERIV
        itrmax = 0
        limit_slope = .false.

        select case (nitr)
        case (0)
           method = NO_DERIV
           itrmax = 0
           limit_slope = .false.
        case (1)
           method = INT_MM
           itrmax = 1
           limit_slope = .false.
        case (2)
           method = INT_MM
           itrmax = 2
           limit_slope = .true.
        case (3)
           method = INT_EMM
           itrmax = 1
           limit_slope = .true.
        case (4)
           method = INT_LV
           itrmax = 1
           limit_slope = .true.
           lv_weight = TWO
        case (5)
           method = INT_NL
           itrmax = 1
           limit_slope = .false.
        case (6,8)                
           method = INT_LV
           itrmax = 1
           limit_slope = .true.
           lv_weight = 1.5_REAL64
        case (7)                  
           method = INT_LV
           itrmax = 1
           limit_slope = .true.
           lv_weight = sqrt(3.0_REAL64)  !JV this value is reduced in inside_com4
        case default
           dimerror = .true.
        end select

        if (dimerror) then
           if (mype.eq.iope) write(*,*)'$$$ DERIVATIVES: nitr = ',nitr
           call global_error('DERIVATIVES: bad value for nitr')
        endif

        if (cell_dim .ne. mesh%cells%numcell_clone) then
           call global_error( &
                & 'derivatives_common_split: cell_dim must be equal to numcell_clone')
        endif

      end subroutine deriv_details
      ! ------------------------------------------------------------------------------
    end subroutine derivatives_common_internal_split
  end module my_derivatives
