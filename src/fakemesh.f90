! Generate a mesh from a PIO file that can be fed to derivatives

module fakemesh
  use mesh_types, only: mesh_t
  use mesh_state_types, only: mesh_state_frac_core_t
  implicit none
  private
  public fakemesh_t
  type :: fakemesh_t
     integer :: ID = 1
     integer :: mpi_id
     integer :: nprocs
     type(mesh_t) :: m
     type(mesh_state_frac_core_t) :: frac_core
   contains
     procedure :: init_PIO_faces
     procedure :: init_PIO_frac_core
     procedure :: init_from_PIO
     procedure :: release_PIO

  end type fakemesh_t

  interface read_and_clone
     procedure read_and_clone_r64
     procedure read_and_clone_i64
     procedure read_and_clone_i32      
  end interface read_and_clone
     
contains
  subroutine read_and_clone_r64(array, name, pioid, iStart, nCount, index)
    use iso_c_binding, only: c_int
    use iso_fortran_env, only: REAL64, INT64
    use clone_lib_module, only: clone_get
    use pio_interface
    implicit none
    real(REAL64), intent(out) :: array(:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: pioid
    integer(INT64), intent(in) :: iStart, nCount
    integer, optional, intent(in) :: index
    real(REAL64), pointer :: tmp(:)
    integer(c_int) :: i

    if ( present(index) ) then
       i = index
    else
       i = 0
    end if

    tmp => pio_get_range_d(pioid, name, i, iStart, nCount)
    array(1:nCount) = tmp
    call pio_release(tmp)
    call clone_get(array)
  end subroutine read_and_clone_r64
  
  subroutine read_and_clone_i64(array, name, pioid, iStart, nCount, index)
    use iso_c_binding, only: c_int, c_int64_t
    use iso_fortran_env, only: INT64
    use pio_interface
    use clone_lib_module, only: clone_get
    implicit none
    integer(INT64), intent(out) :: array(:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: pioid
    integer(INT64), intent(in) :: iStart, nCount
    integer, optional, intent(in) :: index
    integer(c_int64_t), pointer :: tmp(:)
    integer(c_int) :: i

    if ( present(index) ) then
       i = index
    else
       i = 0
    end if

    tmp => pio_get_range_i64(pioid, name, i, iStart, nCount)
    array(1:nCount) = tmp
    call pio_release(tmp)
    call clone_get(array)
  end subroutine read_and_clone_i64
  
  subroutine read_and_clone_i32(array, name, pioid, iStart, nCount, index)
    use iso_c_binding, only: c_int, c_int64_t
    use iso_fortran_env, only: INT64
    use pio_interface
    use clone_lib_module, only: clone_get
    implicit none
    integer, intent(out) :: array(:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: pioid
    integer(INT64), intent(in) :: iStart, nCount
    integer, optional, intent(in) :: index
    integer(c_int64_t), pointer :: tmp(:)
    integer(c_int) :: i

    if ( present(index) ) then
       i = index
    else
       i = 0
    end if

    tmp => pio_get_range_i64(pioid, name, i, iStart, nCount)
    array(1:nCount) = tmp
    call pio_release(tmp)
    call clone_get(array)
  end subroutine read_and_clone_i32
  
  subroutine init_PIO_frac_core(self)
    ! Initialize the frac_core values from file
    use pio_interface
    use define_kind, only: INT64
    implicit none
    class(fakemesh_t) :: self
    integer(INT64) :: iStart, nCount
    
    ASSOCIATE(                        &
         m => self%m,                 &
         pioid => self%ID,            &
         cells => self%m%cells,       &
         faces => self%m%faces,       &
         nprocs => self%nprocs,       &
         mpiid => self%mpi_id,        &
         frac_core => self%frac_core  &
         )

         iStart = cells%cell_address(mpiid)
         nCount = cells%cell_address(mpiid + 1) - iStart
         !loCell =>  pio_get_range_i64(self%id, "cell_index", 2 * idim - 1, iStart, nCount)

         
         
    END ASSOCIATE
      
    
  end subroutine init_PIO_frac_core
  subroutine release_mesh(m)
    use mem_release, only: release
    type(mesh_t), intent(inout) :: m
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
  end subroutine release_mesh

  subroutine release_PIO(self)
    use mesh_types, only: release_mesh
    use pio_interface
    implicit none
    class(fakemesh_t) :: self
    call release_mesh(self%m)
    call pio_release(self%ID)
  end subroutine release_PIO

  subroutine allocate_mesh_scalars(m)
    use iso_c_binding
    use mesh_types, only: nullify_mesh
    implicit none
    type(mesh_t), intent(out) :: m

    allocate( &
         m%cells, &
         m%faces, &
         m%levels, &
         m%neighbors, &
         m%sim, &
         m%amr_vars)

    call nullify_mesh(m)
    
    allocate(m%cells%numcell, m%cells%sum_numcell, m%cells%max_numcell)
    allocate(m%cells%numcell_clone, m%cells%mxcell)

    
    ASSOCIATE(levels => m%levels)
      allocate(levels%numtop, levels%allnumtop)
    END ASSOCIATE

  end subroutine allocate_mesh_scalars

  function gen_partition(ID, ndim, ncell, nprocs, myID, iStart, nCount) result(values)
    ! If nprocs matches number in file, return original partition, otherwise
    ! generate a new partition based on blocks
    use iso_c_binding
    use pio_interface
    implicit None
    integer, intent(in) :: ID, ndim
    integer(c_int64_t), dimension(:), pointer :: values
    integer(c_int64_t), intent(in) :: ncell
    integer, intent(in) :: nprocs, myID
    integer(c_int64_t), intent(out) :: iStart, nCount

    integer :: oldprocs
    integer(c_int64_t), parameter :: one = 1
    integer(c_int64_t) :: nBlocks, i, blockSize
    integer(c_int64_t), dimension(:), pointer :: oldvalues
    real(c_double) :: quantum, next


    allocate(values(0:nprocs))
    oldprocs = pio_length(ID, "global_numcell")
    if (oldprocs == nprocs ) then
       oldvalues => pio_get_i64(ID, "global_numcell", 0)
       values(0) = 1
       do i = 1, nprocs
          values(i) = oldvalues(i) + values(i-1)
       end do
       call pio_release(oldvalues)
    else
       blockSize = 2 ** ndim
       nBlocks = ncell / blockSize

       quantum = real(nBlocks, kind=c_double) / real(nprocs, kind=c_double)

       do i = 0, nprocs-1
          next = quantum * real(i,kind=c_double)
          values(i) = one + blockSize * int(next, kind=c_int64_t)
       end do
       values(nprocs) = ncell + 1
    end if
    iStart = values(myID) 
    nCount = values(myID+1) - values(myID)
    if (myID == 0) then
       write(*,*) 'MPI Partitioning:'
       do i = 0, nprocs-1
          write(*, *) i, values(i), values(i+1), values(i+1) - values(i)
       end do
       write(*,*) '---------------'
    end if

  end function gen_partition

  subroutine init_PIO_faces(self, iStart, nCount, nbrs)
    ! initialize faces from piofile
    ! Missing low side coarse on high boundary faces
    use define_kind, only: LO_SIDE, HI_SIDE
    use mem_release, only: release
    use iso_c_binding
    use clone_lib_module, only: clone_myid
    use pio_interface
    class(fakemesh_t) :: self
    integer(c_int64_t), intent(in) :: iStart, nCount
    integer(c_int64_t), intent(in), dimension(:,:) :: nbrs
    integer(c_int64_t) :: iEnd
    integer(kind=c_int64_t) :: nFace(5), nFaces(5,3)
    integer :: ndim, idim, ilvl, nFaceTypes, iTmp
    integer(c_int64_t) :: iCell, id_lo,id_hi, iFace, offset_now, maxFaces, iTop
    integer(c_int64_t) :: faceIndex(5)
    integer(c_int64_t) :: idxClone
    integer, dimension(:), pointer :: cell_level  
    integer :: idMap(5,3)  ! Maps real ID to face_id array

    ! Two pass face creation - one pass for counting and one for creating

    ASSOCIATE(m => self%m, faces => self%m%faces)
      iEnd = iStart + nCount - 1
      ndim = m%sim%numdim

      cell_level => m%levels%cell_level

      ! Count faces in all directions
      allocate( faces%face_num(ndim) )           
      nFaces = 0
      maxFaces = 0
      META_DIM: do idim = 1, ndim
         META_CELL: do iTop = 1, m%levels%numtop
            iCell = m%levels%ltop(iTop)
            
            id_lo = nbrs(iCell, 2 * idim - 1)
            iFace = get_face_type(iCell, id_lo, LO_SIDE)
            nFaces(iFace, idim) = nFaces(iFace, idim) + 1

            id_hi = nbrs(iCell, 2 * idim)
            iFace = get_face_type(iCell, id_hi, HI_SIDE)
            nFaces(iFace,idim) = nFaces(iFace,idim) + 1
         end do META_CELL
         ! Total number of faces in this direction
         faces%face_num(idim) = sum(nFaces(:,idim))
         maxFaces = max(maxFaces, faces%face_num(idim))
      end do META_DIM
      
      ! Count max number of face types
      idMap = -1
      nFaceTypes = 0
      do iDim = 1, ndim
         iTmp = 0
         do iFace = 1, 5
            if ( nFaces(iFace, iDim) > 0 ) then
               iTmp = iTmp + 1
               idMap(iFace, iDim) = iTmp
            end if
         end do
         if ( iTmp > nFaceTypes ) then
            nFaceTypes = iTmp
         end if
      end do

      ! Allocate space to hold face data
      allocate(  &
           faces%face_id(nFaceTypes,ndim),     &
           faces%face_lo(nFaceTypes, ndim),    &
           faces%face_hi(nFaceTypes, ndim),    &
           faces%face_local(maxFaces, 2, ndim) &
           )

      ! Populate face meta data
      do idim = 1, ndim
         offset_now = 1
         do iFace = 1, 5
            if ( idMap(iFace, iDim) > 0 ) then
               iTmp = idMap(iFace, iDim)
               faces%face_lo(iTmp, iDim) = offset_now
               faces%face_hi(iTmp, iDim) = offset_now + nFaces(iFace, iDim) - 1
               offset_now = offset_now + nFaces(iFace, iDim)
            end if
         end do
      end do

      ! Populate face data
      LOOP_DIM: do idim = 1, ndim
         faceIndex = 0
         LOOP_CELL: do iTop = 1, m%levels%numtop
            iCell = m%levels%ltop(iTop)
            id_lo = nbrs(iCell, 2 * idim - 1)
            iTmp = idMap(get_face_type(iCell, id_lo, LO_SIDE), iDim)
            iFace = faces%face_lo(iTmp, iDim) + faceIndex(iTmp)
            faces%face_local(iFace, LO_SIDE, idim) = id_lo
            faces%face_local(iFace, HI_SIDE, idim) = iCell
            faceIndex(iTmp) = faceIndex(iTmp) + 1

            id_hi = nbrs(iCell, 2 * idim )
            iTmp = idMap(get_face_type(iCell, id_hi, HI_SIDE), iDim)
            iFace = faces%face_lo(iTmp, iDim) + faceIndex(iTmp)
            faces%face_local(iFace, LO_SIDE, idim) = iCell
            faces%face_local(iFace, HI_SIDE, idim) = id_hi
            faceIndex(iTmp) = faceIndex(iTmp) + 1
         end do LOOP_CELL
      end do LOOP_DIM

    END ASSOCIATE
  contains
    integer function get_face_type(the_cell, the_id, the_side)
      use iso_fortran_env, only: INT64
      implicit none
      integer(INT64), intent(in) :: the_cell, the_id
      integer, intent(in) :: the_side
      if ( the_id == the_cell ) then
         if (the_side == LO_SIDE) then
            get_face_type = 1
         else
            get_face_type = 2
         end if
      else if ( cell_level(the_id) < cell_level(the_cell) ) then
         if (the_side == LO_SIDE) then
            get_face_type = 5
         else
            get_face_type = 4
         end if
      else if ( cell_level(the_id) > cell_level(the_cell) ) then
         if (the_side == LO_SIDE) then
            get_face_type = 4
         else
            get_face_type = 5
         end if
      else 
         get_face_type = 3
      end if
    end function get_face_type
    
  end subroutine init_PIO_faces
  
  subroutine init_from_PIO(self, piofile, mpinprocs, mpiid)
    ! Initializes a mesh from a PIO file
    use iso_fortran_env, only: INT64, REAL64
    use clone_lib_module, only: clone_get, clone_base_init, clone_init
    use pio_interface
    implicit none

    class(fakemesh_t) :: self
    character(len=*) :: piofile
    integer, intent(in), optional :: mpinprocs
    integer, intent(in), optional :: mpiid
    integer(c_int64_t), pointer, dimension(:) :: daughter
    integer(c_int64_t) :: i, j, iStart, nCount, myProcs, nCell, iCell
    integer :: nprocs, myid, ndim, iTmp, iDim
    real(c_double), pointer, dimension(:) :: tmp_d
    integer(c_int64_t), dimension(:), pointer :: lo_cell, hi_cell
    integer(c_int64_t), dimension(:), pointer :: amhc_i
    real(c_double), pointer, dimension(:,:) :: cell_hi_lo
    integer(INT64), allocatable, dimension(:,:) :: nbrs


    !*-- Get the MPI numprocs and our processor ID
    
    if (present(mpinprocs)) then
       nprocs = mpinprocs
       if (present(mpiid)) then
          myid = mpiid
       else
          myid = 0
       end if
    else
       call clone_base_init(myid, nprocs)
    end if


    self%mpi_id = myid
    self%nprocs = nprocs
    
    ASSOCIATE( m => self%m, pioid => self%id )

      
      !*-- Initialize the PIO class
      if (myid == 0) write(*,*) 'reading PIO: ', piofile

      !*-- TODO(sriram): Change pio_init() to do a "bare" initialization
      !*--               where it does not read in the cell and daughter
      !*--               arrays for whole simulation
      call pio_init_par(pioid, piofile, nprocs, myid, 0, 1)
      
      !*-- Read in the total number of cells and dimensions
      nCell = pio_ncell(pioid)
      ndim = pio_ndim(pioid)
      
      !*-- Allocate scalars in mesh data structure
      call allocate_mesh_scalars(m)

      !*-- Set ndim
      m%sim%numdim = nDim

      !*-- Generate the MPI partitioning, and current PE's iStart and nCount
      m%cells%cell_address(0:nprocs) => &
           gen_partition(PIOID, ndim, nCell, nprocs, myID, iStart, nCount)

      !*-- Allocate mesh scalars 
      m%cells%numcell = nCount
      m%cells%numcell_clone = nCount
      m%cells%sum_numcell = nCount
      m%cells%max_numcell = nCount

      !*-- Read in neighbors for face and clone processing
      allocate(nbrs(nCount, 2 * ndim))
      do idim = 1, ndim
         lo_Cell => pio_get_range_i64(self%id, "cell_index", 2 * idim - 1, iStart, nCount)
         nbrs(:,2 * idim - 1) = lo_Cell
         call pio_release(lo_Cell)
         nullify(lo_Cell)

         hi_Cell => pio_get_range_i64(self%id, "cell_index", 2 * idim, iStart, nCount)
         nbrs(1:nCount,2 * idim) = hi_Cell
         call pio_release(hi_Cell)
         nullify(hi_Cell)
      end do

      !*-- Count number of top level cells
      daughter => pio_daughter(self%id)
      daughter => pio_get_range_i64(self%id, "cell_daughter", 0, iStart, nCount)
      m%levels%numtop = 0
      m%levels%allnumtop = 0
      do i =1, m%cells%numcell
         if (daughter(i) <= 0) then 
            m%levels%numtop = m%levels%numtop  + 1
         end if
      end do

      !*-- Initialize ltop
      allocate(m%levels%ltop(m%levels%numtop))
      iTmp = 0
      do i = 1, m%cells%numcell
         if (daughter(i) <= 0) then
            iTmp = iTmp + 1
            m%levels%ltop(iTmp) = i
         end if
      end do
      call pio_release(daughter)
      
      !*-- initialize clones
      call clone_init(self%m, nbrs, myid, m%cells%cell_address(0:nprocs))
      
      !*-- Update cell centers
      allocate(m%cells%cell_center(m%cells%numcell_clone, ndim))
      if (myid == 0) write(*,*) 'getting centers'
      do iDim = 1, ndim
         call read_and_clone(m%cells%cell_center(:,iDim), "cell_center", self%id, iStart, nCount, iDim)
      end do
      
      !*-- Update volumes
      allocate(m%cells%vcell(m%cells%numcell_clone))
      call read_and_clone(m%cells%vcell, "vcell", self%id, iStart, nCount)

      !*-- Set high and low half volumes assuming cartesian grid
      allocate(m%cells%cell_half_hi(m%cells%numcell_clone, ndim))
      allocate(m%cells%cell_half_lo(m%cells%numcell_clone, ndim))
      do iDim = 1, ndim
         m%cells%cell_half_lo(:, iDim) = m%cells%vcell/2.0_REAL64
         m%cells%cell_half_hi(:, iDim) = m%cells%vcell/2.0_REAL64
      end do
      
      !*-- Update AMR cell_level
      allocate(m%levels%cell_level(m%cells%numcell_clone))
      call read_and_clone(m%levels%cell_level, "cell_level", self%id, iStart, nCount)

      call self%init_PIO_faces(iStart, nCount, nbrs)
      deallocate(nbrs)
    END ASSOCIATE
  end subroutine init_from_PIO

end module fakemesh

