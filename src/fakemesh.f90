! Generate a mesh from a PIO file that can be fed to derivatives

module fakemesh
  use mesh_types, only: mesh_t
  use mesh_state_types, only: mesh_state_frac_core_t
  use clone_lib_module, only: clone_myid
  implicit none
  public
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
  
  subroutine init_PIO_frac_core(self, iStart, nCount)
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

         ! Get the material counts
         frac_core%vol%obj = pio_get_range_matvar(self%ID, "chunk_vol", 0_INT64, nCount) 

         
    END ASSOCIATE
      
    
  end subroutine init_PIO_frac_core

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
    allocate(m%levels%numtop, m%levels%allnumtop)

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
    ! if (myID == 0) then
    !    write(*,*) 'MPI Partitioning:'
    !    do i = 0, nprocs-1
    !       write(*, *) i, values(i), values(i+1), values(i+1) - values(i)
    !    end do
    !    write(*,*) '---------------'
    ! end if

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
    integer(kind=c_int64_t) :: nFace(5), nFaces(5,3)
    integer :: ndim, idim, ilvl, nFaceTypes, iTmp, iIndex, n_shift, jTmp
    integer(c_int64_t) :: iCell, id_lo,id_hi, iFace, offset_now, maxFaces, iTop, iType
    integer(c_int64_t) :: faceIndex(5)
    integer(c_int64_t) :: idxClone
    integer, dimension(:), pointer :: cell_level  
    integer :: idMap(5,3)  ! Maps real ID to face_id array
    integer, parameter :: offsets_n(3,3) =  reshape([2,4,6, 1,4,5, 1,2,3],[3,3])
    integer :: iCheck


    ! Two pass face creation - one pass for counting and one for creating

    ASSOCIATE(m => self%m, faces => self%m%faces)
      ndim = m%sim%numdim

      if (nDim > 1 ) then
         n_shift = 2 * ndim - 3
      else
         n_shift = 0
      end if
      
      cell_level => m%levels%cell_level

      ! Count faces in all directions
      !
      ! All low side faces belong to the current cell unless
      ! the low side boundary cell is finer (type 4 face).
      !
      ! High side faces count when on physical / PE boundaries
      ! and when on an AMR boundary with a coarse cell on the
      ! high side (type-4 face)

      nFaces = 0
      maxFaces = 0
      META_DIM: do idim = 1, ndim
         META_CELL: do iTop = 1, m%levels%numtop
            iCell = m%levels%ltop(iTop)

            id_lo = nbrs(iCell, 2 * idim - 1)
            iType = get_face_type(iCell, id_lo, LO_SIDE, m%cells%numcell)
            if (iType == 5 .or. iType == 3 .or. iType == 1) then
               nFaces(iType, idim) = nFaces(iType, idim) + 1
               if (iType == 5 .and. id_lo > m%cells%numcell) then
                  ! Add coplanar Off-processor AMR faces
                  nFaces(iType, idim) = nFaces(iType, idim) + n_shift
               end if
            end if
               
            id_hi = nbrs(iCell, 2 * idim)
            iType = get_face_type(iCell, id_hi, HI_SIDE, m%cells%numcell)
            if ((iType == 2 .or. iType == 4) .and. id_hi <= m%cells%numcell) then
               ! Add Face
               nFaces(iType,idim) = nFaces(iType,idim) + 1
            else if (id_hi > m%cells%numcell) then
               ! Need to add PE boundary faces
               nFaces(iType,idim) = nFaces(iType,idim) + 1
            end if
         end do META_CELL
         ! Total number of faces in this direction
         maxFaces = max(maxFaces, sum(nFaces(:,idim)))
      end do META_DIM

      ! Count max number of face types
      idMap = -1
      nFaceTypes = 0
      do iDim = 1, ndim
         iIndex = 0
         do iType = 1, 5
            if ( nFaces(iType, iDim) > 0 ) then
               iIndex = iIndex + 1
               idMap(iType, iDim) = iIndex
            end if
         end do
         if ( iIndex > nFaceTypes ) then
            nFaceTypes = iIndex
         end if
      end do

      ! Allocate space to hold face data
      allocate(  &
           faces%face_num(ndim),               &          
           faces%face_id(nFaceTypes,ndim),     &
           faces%face_lo(nFaceTypes, ndim),    &
           faces%face_hi(nFaceTypes, ndim),    &
           faces%face_local(maxFaces, 2, ndim) &
           )

      ! Populate face meta data
      faces%face_num = 0
      faces%face_local = -1
      faces%face_lo = -1
      faces%face_hi = -1
      do idim = 1, ndim
         offset_now = 1
         do iType = 1, 5
            if ( idMap(iType, iDim) > 0 ) then
               faces%face_num(idim) = faces%face_num(idim) + 1
               iIndex = idMap(iType, iDim)
               faces%face_id(iIndex, iDim) = iType
               faces%face_lo(iIndex, iDim) = offset_now
               faces%face_hi(iIndex, iDim) = offset_now + nFaces(iType, iDim) - 1
               offset_now = offset_now + nFaces(iType, iDim)
            end if
         end do
      end do

      ! Populate face data
      LOOP_DIM: do idim = 1, ndim
         faceIndex = 0
         LOOP_CELL: do iTop = 1, m%levels%numtop
            iCell = m%levels%ltop(iTop)

            ! Low side
            id_lo = nbrs(iCell, 2 * idim - 1)
            iType = get_face_type(iCell, id_lo, LO_SIDE, m%cells%numcell)
            if (iType == 5 .or. iType == 3 .or. iType == 1) then
               iIndex = idMap(iType, iDim)
               iFace = faces%face_lo(iIndex, iDim) + faceIndex(iIndex)
               faces%face_local(iFace, LO_SIDE, idim) = id_lo
               faces%face_local(iFace, HI_SIDE, idim) = iCell
               faceIndex(iIndex) = faceIndex(iIndex) + 1
               if (iType == 5 .and. id_lo > m%cells%numcell) then
                  ! Add in n_shift more faces for off-processor AMR
                  ! do nothing for off processor for now
               end if
            end if

            
            ! High side
            id_hi = nbrs(iCell, 2 * idim )
            iType = get_face_type(iCell, id_hi, HI_SIDE, m%cells%numcell)
            if ((iType == 2 .or. iType == 4) .and. id_hi <= m%cells%numcell) then
               ! Need to add PE boundary faces
               iIndex = idMap(iType, iDim)
               iFace = faces%face_lo(iIndex, iDim) + faceIndex(iIndex)
               faces%face_local(iFace, LO_SIDE, idim) = iCell
               faces%face_local(iFace, HI_SIDE, idim) = id_hi
               faceIndex(iIndex) = faceIndex(iIndex) + 1
            else if (id_hi > m%cells%numcell) then
               ! Need to add PE boundary faces
               iIndex = idMap(iType, iDim)
               iFace = faces%face_lo(iIndex, iDim) + faceIndex(iIndex)
               faces%face_local(iFace, LO_SIDE, idim) = iCell
               faces%face_local(iFace, HI_SIDE, idim) = id_hi
               faceIndex(iIndex) = faceIndex(iIndex) + 1
               ! if (id_hi > m%cells%numcell) then
               !    ! Only triggered on PE boundaries
               !    ! Add in n_shift more faces
               !    do jTmp = 1, n_shift
               !       iFace = faces%face_lo(iIndex, iDim) + faceIndex(iIndex)
               !       faces%face_local(iFace, LO_SIDE, idim) = iCell
               !       faces%face_local(iFace, HI_SIDE, idim) = id_hi + jTmp
               !       faceIndex(iIndex) = faceIndex(iIndex) + 1
               !       iCheck = faces%face_local(iFace, HI_SIDE, idim)
               !    end do
               ! end if
            end if
         end do LOOP_CELL
      end do LOOP_DIM

    END ASSOCIATE
  contains
    pure integer function get_face_type(the_cell, the_nbr, the_side, numcell)
      use iso_fortran_env, only: INT64
      implicit none
      integer(INT64), intent(in) :: the_cell, the_nbr
      integer, intent(in) :: the_side, numcell

      if (the_nbr > numcell) then
         ! Force off-processor to type 3
         get_face_type = 3
      else if ( the_nbr == the_cell ) then
         if (the_side == LO_SIDE) then
            get_face_type = 1
         else
            get_face_type = 2
         end if
      else if ( cell_level(the_nbr) < cell_level(the_cell) ) then
         if (the_side == LO_SIDE) then
            get_face_type = 5
         else
            get_face_type = 4
         end if
      else if ( cell_level(the_nbr) > cell_level(the_cell) ) then
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
    use clone_lib_module, only: mycomm, clone_get, clone_base_init, clone_init, clone_barrier
    use pio_interface
    implicit none

    class(fakemesh_t) :: self
    character(len=*) :: piofile
    integer, intent(in), optional :: mpinprocs
    integer, intent(in), optional :: mpiid
    integer(c_int64_t), pointer, dimension(:) :: daughter
    integer(c_int64_t), pointer, dimension(:) :: mylevel
    integer(c_int64_t) :: i, j, iStart, nCount, myProcs, nCell, iCell, iNbr, myNbr
    integer :: nprocs, myid, ndim, iTmp, iDim
    real(c_double), pointer, dimension(:) :: tmp_d
    integer(c_int64_t), dimension(:), pointer :: lo_cell, hi_cell
    integer(c_int64_t), dimension(:), pointer :: amhc_i
    real(c_double), pointer, dimension(:,:) :: cell_hi_lo
    integer(INT64), allocatable, dimension(:,:) :: nbrs
    integer(INT64), parameter :: offset_n(3) = (/1,2,4/)
    integer(c_int64_t) :: iEnd


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

      !*-- initialize a bare "parallel" PIO struct
      call pio_init_par(pioid, piofile, 0, 1, mycomm)
      
      !*-- Read in the total number of cells and dimensions
      nCell = pio_ncell(pioid)
      ndim = pio_ndim(pioid)
      
      !*-- Allocate scalars in mesh data structure
      call allocate_mesh_scalars(m)

      !*-- Set ndim
      m%sim%numdim = nDim

      !*-- Generate the MPI partitioning, and current PE's iStart and nCount
      m%cells%cell_address => &
           gen_partition(PIOID, ndim, nCell, nprocs, myID, iStart, nCount)
      iEnd = iStart + nCount -1

      !SS call pio_init_materials(pioid, iStart, nCount)
      
      !*-- Allocate mesh scalars 
      m%cells%numcell = nCount
      m%cells%numcell_clone = nCount
      m%cells%sum_numcell = nCount
      m%cells%max_numcell = nCount

      !*-- Read in neighbors for face and clone processing
      allocate(nbrs(nCount, 2 * ndim))
      nbrs = -20
      do idim = 1, ndim
         lo_Cell => pio_get_range_i64(self%id, "cell_index", 2 * idim - 1, iStart, nCount)
         nbrs(1:nCount,2 * idim - 1) = lo_Cell(1:nCount)
         call pio_release(lo_Cell)
         nullify(lo_Cell)

         hi_Cell => pio_get_range_i64(self%id, "cell_index", 2 * idim, iStart, nCount)
         nbrs(1:nCount,2 * idim) = hi_Cell(1:nCount)
         call pio_release(hi_Cell)
         nullify(hi_Cell)
      end do

      iEnd = iStart + nCount - 1
      !*-- Count number of top level cells
      ! Fix neighbors array for refined neighbors
      m%levels%cell_daughter => pio_get_range_i64(self%id, "cell_daughter", 0, iStart, nCount)
      daughter => m%levels%cell_daughter
      m%levels%numtop = 0
      m%levels%allnumtop = 0
      do i =1, m%cells%numcell
         if (daughter(i) <= 0) then 
            m%levels%numtop = m%levels%numtop  + 1
            do iDim = 1, nDim
               ! Check low side for type 4 face
               iNbr = nbrs(i, 2*iDim-1)
               myNbr = iNbr - iStart + 1
               if (iNbr >= iStart .and. iNbr <= iEnd) then
                  if (daughter(myNbr) > 0 ) then
                     nbrs(i,2*iDim-1) = daughter(myNbr) + offset_n(iDim)
                  end if
               end if
               
               ! Check high side for type 5 face
               iNbr = nbrs(i, 2*iDim)
               myNbr = iNbr - iStart + 1
               if (iNbr >= iStart .and. iNbr <= iEnd) then
                  if (daughter(myNbr) > 0) then
                     nbrs(i,2*iDim) = daughter(myNbr)
                  end if
               end if
            end do
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

      !*-- initialize clones
      call clone_init(self%m, nbrs, pioID)

      call clone_barrier()
      if (myid == 0 ) write(*,*) 'Done initializing, reading data'
      !*-- Update cell centers
      allocate(m%cells%cell_center(m%cells%numcell_clone, ndim))
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

      !*-- update the mesh state variables
      !SS call self%init_PIO_frac_core(iStart, nCount)
      
      if (myid == 0 ) write(*,*) 'Done reading data, initializing faces'
      call self%init_PIO_faces(iStart, nCount, nbrs)
      deallocate(nbrs)

      call clone_barrier()
      if (myid == 0 ) write(*,*) 'Done initializing faces'
    END ASSOCIATE
  end subroutine init_from_PIO

end module fakemesh

