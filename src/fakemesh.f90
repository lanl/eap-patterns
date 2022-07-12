! Generate a mesh from a PIO file that can be fed to derivatives

module fakemesh
  use mesh_types, only: mesh_t
  implicit none
  private
  public fakemesh_t
  type :: fakemesh_t
     integer :: ID = 1
     type(mesh_t) :: m
   contains
     procedure :: init_PIO_faces
     procedure :: init_from_PIO
     procedure :: release_PIO

  end type fakemesh_t
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

  subroutine allocate_base_mesh(m)
    use iso_c_binding
    implicit none
    type(mesh_t), intent(out) :: m

    allocate( &
         m%cells, &
         m%faces, &
         m%levels, &
         m%neighbors, &
         m%sim, &
         m%amr_vars)

    allocate(m%cells%numcell, m%cells%sum_numcell, m%cells%max_numcell)
    allocate(m%cells%numcell_clone, m%cells%mxcell)

    ASSOCIATE(levels => m%levels)
      allocate(levels%numtop, levels%allnumtop)
    END ASSOCIATE

  end subroutine allocate_base_mesh

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
    nCount = values(myID+1) - values(myID) - 1
    write(*,*) 'MPI Partitioning:'
    do i = 1, nprocs
       write(*, *) i, values(i-1), values(i), values(i) - values(i-1)
    end do
    write(*,*) '---------------'

  end function gen_partition

  subroutine init_PIO_faces(self, iStart, nCount)
    ! initialize faces from piofile
    use define_kind, only: LO_SIDE, HI_SIDE
    use iso_c_binding
    use pio_interface
    class(fakemesh_t) :: self
    integer(c_int64_t) :: iStart, nCount
    integer(kind=c_int64_t) :: nFace(5), nFaces(5,3)
    integer(c_int64_t), dimension(:), pointer :: cell_level     ! The levels of the cells, i32 should be fine...
    integer(c_int64_t), dimension(:), pointer :: loCell, hiCell ! The low and high side cells
    integer :: ndim, idim, ilvl, nFaceTypes, iTmp
    integer(c_int64_t) :: iCell, id_lo,id_hi, iFace, offset_now, maxFaces
    integer(c_int64_t), allocatable, dimension(:,:) :: nbrs
    integer(c_int64_t) :: faceIndex(5)
    integer :: idMap(5,3)  ! Maps real ID to face_id array

    ! Two pass face creation - one pass for counting and one for creating

    write(*,*) 'initializing faces'
    ASSOCIATE(m => self%m, faces => self%m%faces)
      ndim = pio_ndim(self%id)


      write(*,*) 'ndim=', ndim, 'ncount=', nCount
      allocate(nbrs(nCount, 2 * ndim))
      ! reading whole array to avoind having to figure out
      ! which cells to read
      cell_level => pio_get_i64(self%id, "cell_level", 0) 


      ! Count faces in all directions
      allocate( faces%face_num(ndim) )           
      nFaces = 0
      maxFaces = 0
      META_DIM: do idim = 1, ndim
         loCell =>  pio_get_range_i64(self%id, "cell_index", 2 * idim - 1, iStart, nCount) 
         nbrs(:,2 * idim - 1) = loCell(:nCount)
         call pio_release(loCell)
         nullify(loCell)

         hiCell =>  pio_get_range_i64(self%id, "cell_index", 2 * idim, iStart, nCount)
         nbrs(:,2 * idim) = hiCell(:nCount)
         call pio_release(hiCell)
         nullify(hiCell)

         META_CELL: do iCell = 1, nCount
            !write(*,*) iCell, idim
            id_lo = nbrs(iCell, 2 * idim - 1)
            if ( id_lo == iCell + iStart - 1 ) then
               nFaces(1,idim) = nFaces(1,idim) + 1
            else if ( cell_level(id_lo) < cell_level(iCell) ) then
               nFaces(5,idim) = nFaces(5,idim) + 1
            else if ( cell_level(id_lo) > cell_level(iCell) ) then
               nFaces(4,idim) = nFaces(4,idim) + 1
            else 
               nFaces(3,idim) = nFaces(3,idim) + 1
            end if

            id_hi = nbrs(iCell, 2 * idim )
            if ( id_hi == iCell + iStart - 1_c_int64_t ) then
               nFaces(2,idim) = nFaces(1,idim) + 1
            else if ( cell_level(id_hi) < cell_level(iCell) ) then
               nFaces(4,idim) = nFaces(4,idim) + 1
            else if ( cell_level(id_hi) > cell_level(iCell) ) then
               nFaces(5,idim) = nFaces(5,idim) + 1
            else 
               nFaces(3,idim) = nFaces(3,idim) + 1
            end if
         end do META_CELL

         ! Total number of faces in this direction
         faces%face_num(idim) = sum(nFaces(:,idim))
         maxFaces = max(maxFaces, faces%face_num(idim))
      end do META_DIM

      ! Count max number of face types
      idMap = 0
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
         iTmp = 1
         offset_now = 1
         do iFace = 1, 5
            if ( nFaces(iFace, iDim) > 0 ) then
               faces%face_lo(iTmp, iDim) = offset_now
               faces%face_hi(iTmp, iDim) = offset_now + nFaces(iFace, iDim) - 1
               iTmp = iTmp + 1
               offset_now = offset_now + nFaces(iFace, iDim)
            end if
         end do
      end do

      ! Populate face data
      LOOP_DIM: do idim = 1, ndim
         faceIndex = 0
         LOOP_CELL: do iCell = 1, nCount
            !write(*,*) iCell, idim
            id_lo = nbrs(iCell, 2 * idim - 1)
            if ( id_lo == iCell + iStart - 1 ) then
               iTmp = idMap(1,idim)
            else if ( cell_level(id_lo) < cell_level(iCell) ) then
               iTmp = idMap(5,idim)
            else if ( cell_level(id_lo) > cell_level(iCell) ) then
               iTmp = idMap(1,idim)
            else 
               iTmp = idMap(1,idim)
            end if
            iFace = faces%face_lo(iTmp, iDim) + faceIndex(iTmp)
            faces%face_local(iFace, LO_SIDE, idim) = id_lo
            faces%face_local(iFace, HI_SIDE, idim) = iCell
            faceIndex(iTmp) = faceIndex(iTmp) + 1

            id_hi = nbrs(iCell, 2 * idim )
            if ( id_hi == iCell + iStart - 1_c_int64_t ) then
               iTmp = idMap(2,idim)
            else if ( cell_level(id_hi) < cell_level(iCell) ) then
               iTmp = idMap(4,idim)
            else if ( cell_level(id_hi) > cell_level(iCell) ) then
               iTmp = idMap(5,idim)
            else 
               iTmp = idMap(3,idim)
            end if
            iFace = faces%face_lo(iTmp, iDim) + faceIndex(iTmp)
            faces%face_local(iFace, LO_SIDE, idim) = iCell
            faces%face_local(iFace, HI_SIDE, idim) = id_hi
            faceIndex(iTmp) = faceIndex(iTmp) + 1
         end do LOOP_CELL
      end do LOOP_DIM

      call pio_release(cell_level)
      deallocate(nbrs)
    END ASSOCIATE
  end subroutine init_PIO_faces

  subroutine init_from_PIO(self, piofile, mpinprocs, mpiid)
    ! Initializes a mesh from a PIO file
    use pio_interface
    implicit none

    class(fakemesh_t) :: self
    character(len=*) :: piofile
    integer, intent(in), optional :: mpinprocs
    integer, intent(in), optional :: mpiid
    integer(c_int64_t), pointer, dimension(:) :: daughter
    integer(c_int64_t) :: i, j, iStart, nCount, myProcs, nCell, iCell
    integer :: nprocs, myid, ndim
    real(c_double), pointer, dimension(:) :: tmp_d
    integer(c_int64_t), dimension(:), pointer :: cell_level
    integer(c_int64_t), dimension(:), pointer :: amhc_i
    real(c_double), pointer, dimension(:,:) :: cell_hi_lo
    integer :: numlev

    if (present(mpinprocs)) then
       nprocs = mpinprocs
    else
       nprocs = 1
    end if

    if (present(mpiid)) then
       myid = mpiid
    else
       myid = 0
    end if

    ASSOCIATE( m => self%m, id => self%id )
      write(*,*) 'reading PIO: ', piofile
      call pio_init(id, piofile, 1)
      iStart = 1
      nCell = pio_ncell(id)
      ndim = pio_ndim(id)
      write(*,*) 'ncells = ', nCell
      write(*,*) 'nmat = ', pio_nmat(id)

      call allocate_base_mesh(m)
      ! Read in numdim
      m%sim%numdim = pio_ndim(id)

      ! Set up cells
      m%cells%cell_address(0:nprocs) => &
           gen_partition(ID, ndim, nCell, nprocs, myID, iStart, nCount)
      m%cells%numcell_clone = pio_ncell(id)
      m%cells%numcell = m%cells%numcell_clone
      m%cells%sum_numcell = m%cells%numcell_clone
      m%cells%max_numcell = m%cells%numcell_clone

      ! Read in cell centers
      allocate(m%cells%cell_center(nCount, ndim))
      do i = 1, ndim
         tmp_d => pio_get_range_d(id, "cell_center", int(i,kind=c_int), iStart, nCount)
         m%cells%cell_center(1:nCount,i) = tmp_d(1:nCount)
         call pio_release(tmp_d)
      end do

      ! Read in volumes and fill in cell_half_lo/hi
      m%cells%vcell => pio_get_range_d(id, "vcell", 0, iStart, nCount)

      ! Generate cell sizes by level
      cell_level => pio_get_range_i64(id, "cell_level", 0, iStart, nCount)
      amhc_i => pio_get_range_i64(id, "amhc_i", 0, iStart, nCount)
      numlev = amhc_i(45)
      if ( numlev < maxval(cell_level)) then
         numlev = maxval(cell_level)
      end if
      write(*,*) 'numlev = ', numlev
      allocate(cell_hi_lo(numlev, ndim))
      ! We know the size of first cell based on block size
      do i = 1, ndim
         cell_hi_lo(1,i) = m%cells%cell_center(2**ndim-1,i) - m%cells%cell_center(1,i)
      end do
      do i = 2, numlev
         cell_hi_lo(i,:) = cell_hi_lo(i-1,:)/2.0D0
      end do
      allocate(m%cells%cell_half_hi(nCount, ndim), m%cells%cell_half_lo(nCount, ndim))
      do i=1, ndim
         m%cells%cell_half_hi(:,ndim) = cell_hi_lo(cell_level,ndim) / 2.0
         m%cells%cell_half_lo(:,ndim) = cell_hi_lo(cell_level,ndim) /2.0
      end do

      do i = 1, ndim
         write(*,*) 'cell_hilo =', i, cell_hi_lo(:, ndim)
      end do

      ! Set up levels
      daughter => pio_daughter(id)
      m%levels%numtop = 0
      m%levels%allnumtop = 0
      do i =1, m%cells%numcell
         if (daughter(i) <= 0) then 
            m%levels%numtop = m%levels%numtop  + 1
         end if
      end do
      m%levels%allnumtop = m%levels%numtop
      do i =m%cells%numcell+1, m%cells%numcell_clone
         if (daughter(i) <= 0) then 
            m%levels%allnumtop = m%levels%allnumtop  + 1
         end if
      end do

      allocate(&
           m%levels%ltop_nv(m%levels%numtop), &
           m%levels%alltop(m%levels%allnumtop) &
           )
      j = 0
      do i =1, m%cells%numcell
         if (daughter(i) <= 0) then
            j = j + 1
            if (i <= m%cells%numcell) then
               m%levels%ltop_nv(j) = i
            end if
            m%levels%alltop(j) = i
         end if
      end do
      write(*,*) 'cells=',m%cells%numcell,m%cells%numcell_clone,&
           m%levels%numtop, m%levels%allnumtop

      call self%init_PIO_faces(iStart, nCount)
    END ASSOCIATE
  end subroutine init_from_PIO

end module fakemesh

