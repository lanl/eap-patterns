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

! Generate a mesh from a PIO file that can be fed to derivatives

module fakemesh
  use mesh_types, only: mesh_t
  use mesh_state_types, only: mesh_state_frac_core_t, mesh_state_core_t
  use clone_lib_module, only: clone_myid, clone_reduce, CLONE_SUM, clone_abort
  use binreader
  
  implicit none
  public
  public fakemesh_t
  type :: fakemesh_t
     type(binfile) :: bfp
     integer :: mpi_id
     integer :: nprocs
     type(mesh_t) :: m
     type(mesh_state_core_t) :: core
     type(mesh_state_frac_core_t) :: frac_core
   contains
     procedure :: init_faces
     procedure :: init_frac_core
     procedure :: init
     procedure :: release

  end type fakemesh_t

  integer, private :: ncell = -1
contains
  
  subroutine init_frac_core(self, iStart, nCount)
    ! Initialize the frac_core values from file
    use iso_fortran_env, only : INT64, REAL64
    use matdefcm, only: nummat
    implicit none
    class(fakemesh_t) :: self
    integer(INT64) :: iStart, nCount
    real(REAL64), dimension(:,:), pointer :: dptr
    if (clone_myid() == 0 ) write(*,*) '   Initializing_frac_core'
    
    ASSOCIATE(                        &
         m => self%m,                 &
         bfp => self%bfp,            &
         cells => self%m%cells,       &
         faces => self%m%faces,       &
         nprocs => self%nprocs,       &
         mpiid => self%mpi_id,        &
         frac_core => self%frac_core  &
         )

      self%frac_core%vol%nmat = nummat
      self%frac_core%vol%ncells = m%cells%numcell_clone
      allocate(self%frac_core%vol%obj(m%cells%numcell_clone, nummat))
      dptr => self%frac_core%vol%obj
      call self%bfp%read(dptr, "frac_vol", iStart, nCount, int(nummat, kind=INT64))

      self%frac_core%eng%nmat = nummat
      self%frac_core%eng%ncells = m%cells%numcell_clone
      allocate(self%frac_core%eng%obj(m%cells%numcell_clone, nummat))
      dptr => self%frac_core%eng%obj
      call self%bfp%read(dptr, "frac_eng", iStart, nCount, int(nummat, kind=INT64))

      self%frac_core%mass%nmat = nummat
      self%frac_core%mass%ncells = m%cells%numcell_clone
      allocate(self%frac_core%mass%obj(m%cells%numcell_clone, nummat))
      dptr => self%frac_core%mass%obj
      call self%bfp%read(dptr, "frac_mass", iStart, nCount, int(nummat, kind=INT64))


    END ASSOCIATE
      
    
  end subroutine init_frac_core

  subroutine release(self)
    use mesh_types, only: release_mesh
    implicit none
    class(fakemesh_t) :: self
    call release_mesh(self%m)
    call self%bfp%release()
  end subroutine release

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

  function gen_partition(bfp, nprocs, myID, iStart, nCount) result(values)
    ! If nprocs matches number in file, return original partition, otherwise
    ! generate a new partition based on blocks
    use iso_c_binding
    implicit None
    type(binFile) :: bfp
    integer(c_int64_t), dimension(:), pointer :: values
    integer, intent(in) :: nprocs, myID
    integer(c_int64_t), intent(out) :: iStart, nCount

    integer(c_int64_t), parameter :: one = 1
    integer(c_int64_t) :: nBlocks, i, blockSize
    real(c_double) :: quantum, next
    integer(c_INT64_t) :: ndim

    ndim = bfp%ndim
    ncell = bfp%ncells


    allocate(values(0:nprocs))
    blockSize = 2 ** ndim
    nBlocks = ncell / blockSize
    
    quantum = real(nBlocks, kind=c_double) / real(nprocs, kind=c_double)
    
    do i = 0, nprocs-1
       next = quantum * real(i,kind=c_double)
       values(i) = one + blockSize * int(next, kind=c_int64_t)
    end do
    values(nprocs) = ncell + 1

    iStart = values(myID) 
    nCount = values(myID+1) - values(myID)

  end function gen_partition

  subroutine init_faces(self, iStart, nCount, nbrs, face_type, nFaces)
    ! initialize faces
    ! Missing low side coarse on high boundary faces
    use iso_fortran_env, only: INT8, INT64
    use define_kind, only: LO_SIDE, HI_SIDE
    use mem_release, only: release
    use iso_c_binding
    class(fakemesh_t) :: self
    integer(INT64), intent(in) :: iStart, nCount
    integer(INT64), intent(in), dimension(:,:) :: nbrs
    integer(INT8), intent(in) :: face_type(:,:)
    integer(INT64), intent(in) :: nFaces(5,3)
    integer(kind=INT64) :: nFace(5)
    integer :: ndim, idim, ilvl, nFaceTypes, iTmp, iIndex, n_shift, jTmp, j, idx
    integer(INT64) :: iCell, id_lo,id_hi, iFace, offset_now, maxFaces, iTop, iType
    integer(INT64) :: faceIndex(5)
    integer(INT64) :: idxClone
    integer, dimension(:), pointer :: cell_level  
    integer :: idMap(5,3)  ! Maps real ID to face_id array
    integer, parameter :: offsets_n(3,3) =  reshape([2,4,6, 1,4,5, 1,2,3],[3,3])
    integer(INT64) :: face_index(5,3)


    ! Two pass face creation - one pass for counting and one for creating

    ASSOCIATE(m => self%m, faces => self%m%faces)
      ndim = m%sim%numdim

      if (nDim > 1 ) then
         n_shift = 2 * ndim - 3
      else
         n_shift = 0
      end if
      
      cell_level => m%levels%cell_level
    
      ! Count max number of face types
      idMap = -1
      face_index = 0
      maxFaces = 0
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
         maxFaces = max(maxFaces, sum(nFaces(:,idim)))
      end do

      ! Allocate space to hold face data
      allocate(  &
           faces%face_flag(5,3),               &          
           faces%face_num(ndim),               &          
           faces%face_id(nFaceTypes,ndim),     &
           faces%face_lo(nFaceTypes, ndim),    &
           faces%face_hi(nFaceTypes, ndim),    &
           faces%face_local(maxFaces, 2, ndim) &
           )

      ! Populate face meta data
      faces%face_flag = 0
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
               offset_now = offset_now + nFaces(iType, iDim)
               faces%face_hi(iIndex, iDim) = offset_now -1
            end if
         end do
      end do

      ! Populate face data
      LOOP_DIM: do idim = 1, ndim
         faceIndex = 0
         LOOP_CELL: do iTop = 1, m%levels%numtop
            iCell = m%levels%ltop(iTop)
            ! Low side
            idx = 2 * idim - 1
            id_lo = nbrs(idx, iCell)
            iType = face_type(idx, iCell)
            if (iType /= 4 .or. id_lo > m%cells%numcell) then
               iIndex = idMap(iType, iDim)
               iFace = faces%face_lo(iIndex, iDim) + faceIndex(iIndex)
               faces%face_local(iFace, LO_SIDE, idim) = id_lo
               faces%face_local(iFace, HI_SIDE, idim) = iCell
               faceIndex(iIndex) = faceIndex(iIndex) + 1
               if (iType == 4) then
                  ! Only true if id_lo > numcell
                  ! Add in n_shift more faces for off-processor AMR
                  do j = 1, n_shift
                     iFace = faces%face_lo(iIndex, iDim) + faceIndex(iIndex)
                     faces%face_local(iFace, LO_SIDE, idim) = id_lo + 1
                     faces%face_local(iFace, HI_SIDE, idim) = iCell
                     faceIndex(iIndex) = faceIndex(iIndex) + 1
                  end do
               end if
            end if

            
            ! High side
            idx = 2 * idim
            id_hi = nbrs(idx, iCell)
            iType = face_type(idx, iCell)
            iIndex = idMap(iType, iDim)
            iFace = faces%face_lo(iIndex, iDim) + faceIndex(iIndex)
22          FORMAT('icell=',i8,', itype=',i2, ', idx=', i2, ', dim=', i1, ', map=', 5(i2,','))
            if (iType == 2 .or. iType == 4 .or. id_hi > m%cells%numcell) then
               iFace = faces%face_lo(iIndex, iDim) + faceIndex(iIndex)
               faces%face_local(iFace, LO_SIDE, idim) = iCell
               faces%face_local(iFace, HI_SIDE, idim) = id_hi
               faceIndex(iIndex) = faceIndex(iIndex) + 1
               if (iType == 5) then
                  ! Need to add shifted faces
                  ! Only true if id_hi > numcell
                  do jTmp = 1, n_shift
                     iFace = faces%face_lo(iIndex, iDim) + faceIndex(iIndex)
                     faces%face_local(iFace, LO_SIDE, idim) = iCell
                     faces%face_local(iFace, HI_SIDE, idim) = id_hi + 1
                     faceIndex(iIndex) = faceIndex(iIndex) + 1
                  end do
               end if
            end if
         end do LOOP_CELL
      end do LOOP_DIM

    END ASSOCIATE
    
  end subroutine init_faces
  
  subroutine findCloneID(idNew, iNbr, num_clone, clone_map, iType, idx, n_shift, idim)
    use iso_fortran_env, only: INT64
    implicit none
    integer(INT64), intent(out) :: idNew
    integer(INT64), intent(in) :: iNbr
    integer(INT64), intent(inout) :: num_clone, clone_map(:)
    integer, intent(in) :: iType, idx
    integer(INT64), intent(in) :: n_shift
    integer, intent(in) :: idim
    integer, parameter :: offsets_n(3,3) =  reshape([2,4,6, 1,4,5, 1,2,3],[3,3])

    integer :: j

    idNew = -1

    do j = 1, num_clone
       if (clone_map(j) == iNbr) then
          idNew = j
          exit
       end if
    end do

    if (idNew < 0) then
       num_clone = num_clone + 1
       clone_map(num_clone) = iNbr
       idNew = num_clone
       if ( (iType == 4 .and. mod(idx,2) == 1) .or. &
            (iType == 5 .and. mod(idx,2) == 0)) then
          do j = 1, n_shift
             num_clone = num_clone + 1
             clone_map(num_clone) = iNbr + offsets_n(j, idim)
          end do
       end if
    end if
  end subroutine findCloneID

  subroutine init(self, myfile, mpinprocs, mpiid)
    use iso_fortran_env, only: INT64, REAL64, INT8
    use iso_c_binding
    use matdefcm, only: nummat
    use clone_lib_module, only: mycomm, clone_abort, clone_get, clone_base_init, clone_init, clone_barrier, clone_myid
    implicit none

    class(fakemesh_t) :: self
    character(len=*) :: myfile
    integer, intent(in), optional :: mpinprocs
    integer, intent(in), optional :: mpiid
    integer(INT64), pointer, dimension(:) :: daughter
    integer(INT64), pointer, dimension(:) :: clone_map, tmp_clone_map
    integer(INT64) :: i, j, iStart, nCount, myProcs, nCell, iNbr, myNbr, nClone, iNew, tmp_i64
    integer :: nprocs, myid, ndim, iTmp, iDim, idx, iType
    real(c_double), pointer, dimension(:) :: tmp_d
    integer(INT64), dimension(:), pointer :: lo_cell, hi_cell
    integer(INT64), dimension(:), pointer :: amhc_i
    real(c_double), pointer, dimension(:,:) :: cell_hi_lo
    integer(INT8), pointer, dimension(:,:) :: face_type => NULL()
    integer(INT64), pointer, dimension(:,:) :: nbrs => NULL()
    integer(INT64), parameter :: offset_n(3) = (/1,2,4/)
    integer(c_int64_t) :: iEnd, iCell
    character(len=128) :: tmpChar1
    integer(INT64) :: n_shift
    integer, parameter :: offsets_n(3,3) =  reshape([2,4,6, 1,4,5, 1,2,3],[3,3])
    integer(INT64) :: face_count(5,3)


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
    
    ASSOCIATE( m => self%m )
      
      call self%bfp%init(myfile)

      
      !*-- Read in the total number of cells and dimensions
      ndim = self%bfp%nDim
      if ( nDim > 1) then
         n_shift = 2 * nDim - 3
      else
         n_shift = 0
      end if

      nCell = self%bfp%nCells
      nummat = self%bfp%nMat
      
      
      !*-- Allocate scalars in mesh data structure
      call allocate_mesh_scalars(m)

      !*-- Set ndim
      m%sim%numdim = nDim
      m%sim%numvel = nDim

      !*-- Generate the MPI partitioning, and current PE's iStart and nCount
      m%cells%cell_address => &
           gen_partition(self%bfp, nprocs, myID, iStart, nCount)
      iEnd = iStart + nCount -1

      !SS call pio_init_materials(pioid, iStart, nCount)
      
      !*-- Allocate mesh scalars 
      m%cells%numcell = nCount
      m%cells%numcell_clone = nCount
      m%cells%sum_numcell = nCount
      m%cells%max_numcell = nCount

      !*-- Read in neighbors for face and clone processing

      iEnd = iStart + nCount - 1
      !*-- Count number of top level cells and number of clones
      if (myid == 0 ) then
         write(*,*)  '  Reading basic mesh data'
      end if
      do i=0,nprocs-1
         if (i == myID) then
            ! write(*, '("    Reading on processor: ",i6)') i
            call self%bfp%read(nbrs, "cell_index", iStart, nCount, 2_INT64 * nDim)
            call self%bfp%read(face_type, "face_type", iStart, nCount, 2_INT64 * nDim)
            call self%bfp%read(m%levels%cell_daughter, "cell_daughter", iStart, nCount)
         end if
         call clone_barrier()
      end do
      daughter => m%levels%cell_daughter
      

      ! Find top level cells
      !
      ! Count faces in all directions:
      !   All low side faces belong to the current cell unless
      !   the low side boundary cell is finer (type 4 face).
      !
      !   High side faces count when on physical / PE boundaries
      !   and when on an AMR boundary with a coarse cell on the
      !   high side (type-4 face)
      if (myID == 0 ) write(*, *) '  Counting clones'
      nClone = 0
      face_count = 0
      m%levels%numtop = 0
      m%levels%allnumtop = 0
      do iCell = 1, m%cells%numcell
         if (daughter(iCell) == 0) then
            m%levels%numtop = m%levels%numtop  + 1
            ! Check neighbors for top level cells only
            do iDim = 1, nDim
               ! Check low side face
               idx = 2 * iDim - 1
               iNbr = nbrs(idx, iCell)
               iType = face_type(idx, iCell)
               if (iType /= 4) then
                  ! Add all non-type 4 faces
                  face_count(iType, iDim) = face_count(iType, iDim) + 1
                  if (iNbr < iStart .or. iNbr > iEnd) then
                     nClone = nClone + 1
                  end if
               else if (iNbr < iStart .or. iNbr > iEnd) then
                  ! Off processor, only type 4 faces here
                  ! Type 4 faces require additional clones
                  face_count(iType, iDim) = face_count(iType, iDim) + 1 + n_shift
                  nClone = nClone + 1 + n_shift
               end if

               ! Check high side
               idx = 2 * idim
               iNbr = nbrs(idx, iCell)
               iType = face_type(idx, iCell)
               if (iNbr < iStart .or. iNbr > iEnd) then
                  ! Off processor, *always* add
                  face_count(iType, iDim) = face_count(iType, iDim) + 1
                  nClone = nClone + 1

                  ! Type 5 faces require additional clones
                  if (iType == 5) then
                     nClone = nClone + n_shift
                     face_count(iType, iDim) = face_count(iType, iDim) + n_shift
                  end if
               else if (iType == 2 .or. iType == 4) then
                  ! Add on PE neighbor if it is of type 2 or 4 
                  face_count(iType, iDim) = face_count(iType, iDim) + 1
               end if
            end do
         end if
      end do

      if (m%levels%numtop == 0) then
         call clone_abort('No top level cells.  Reduce number of processors')
      end if

      call clone_barrier()
      if (myID == 0 ) write(*, *) '  Adjusting neighbors'
      !*-- Initialize ltop
      !*-- remap neighbors and configure clone_map
      allocate(tmp_clone_map(nClone))
      tmp_clone_map = -1
      allocate(m%levels%ltop(m%levels%numtop))
      m%cells%numcell_clone = m%cells%numcell + nClone
      iTmp = 0
      nClone = 0
      do iCell = 1, m%cells%numcell
         if (daughter(iCell) <= 0) then
            iTmp = iTmp + 1
            m%levels%ltop(iTmp) = iCell
            ! fix neighbors array
            LOOP_2_NDIM: do idx = 1, 2 * nDim
               ! Check low side for type 4 face
               idim = (idx + 1) / 2
               iNbr = nbrs(idx, iCell)
               iType = face_type(idx, iCell)
               if (iNbr >= iStart .and. iNbr <= iEnd) then
                  nbrs(idx, iCell) = iNbr - iStart + 1
               else
                  ! Add a clone for *all* off-PE cells
                  call findCloneID(iNew, iNbr, nClone, tmp_clone_map, iType, idx, n_shift, idim)
                  nbrs(idx, iCell) = iNew + m%cells%numcell
                  ! Note that above call adds the additional clones for sister cells
               end if
            end do LOOP_2_NDIM
         end if
      end do
      call clone_barrier()
      
      if (myID == 0 ) write(*, *) '  Generating clone map'

      ! Set clone_map to correct size and copy data
      allocate(clone_map(nClone))
      clone_map(1:nClone) = tmp_clone_map(1:nClone)
      deallocate(tmp_clone_map)


      m%cells%numcell_clone = m%cells%numcell + nClone
      m%levels%allnumtop = nClone + m%levels%numtop
      allocate(m%levels%alltop(m%levels%allnumtop))
      m%levels%alltop(1:m%levels%numtop) = m%levels%ltop(1:m%levels%numtop)
      do iCell = m%levels%numtop + 1, m%levels%allnumtop
         m%levels%alltop(1:m%levels%numtop) = iCell
      end do

      ! Create faces
      if (myid == 0 ) write(*,*) '  Initializing faces'
      call self%init_faces(iStart, nCount, nbrs, face_type, face_count)
      
      !*-- initialize clones
      if (myid == 0 ) write(*,*) '  Initializing communications'
      call clone_init(self%m, nbrs, clone_map)

      deallocate(nbrs)
      deallocate(clone_map)
      deallocate(m%levels%cell_daughter)

      call clone_barrier()
      if (myid == 0 ) write(*,*) '  Reading cell info'

      !*-- Update cell daughters
      allocate(m%levels%cell_daughter(m%cells%numcell_clone))
      call self%bfp%read(m%levels%cell_daughter, "cell_daughter", iStart, nCount)
      call clone_get(m%levels%cell_daughter)

      !*-- Update cell levels
      allocate(m%levels%cell_level(m%cells%numcell_clone))
      call self%bfp%read(m%levels%cell_level, "cell_level", iStart, nCount)
      call clone_get(m%levels%cell_level)

      !*-- Update cell centers
      allocate(m%cells%cell_center(m%cells%numcell_clone, ndim))
      call self%bfp%read(m%cells%cell_center, "cell_center", iStart, nCount, int(nDim, kind=INT64))
      do iDim=1, nDim
         call clone_get(m%cells%cell_center(:,iDim))
      end do

      !*-- Update volumes
      allocate(m%cells%vcell(m%cells%numcell_clone))
      call self%bfp%read(m%cells%vcell, "vcelll", iStart, nCount)
      call clone_get(m%cells%vcell)
      !*-- Set high and low half volumes assuming cartesian grid
      allocate(m%cells%cell_half_hi(m%cells%numcell_clone, ndim))
      allocate(m%cells%cell_half_lo(m%cells%numcell_clone, ndim))
      do iDim = 1, ndim
         m%cells%cell_half_lo(:, iDim) = m%cells%vcell/2.0_REAL64
         m%cells%cell_half_hi(:, iDim) = m%cells%vcell/2.0_REAL64
      end do
      !*-- update the mesh state variables
      call self%init_frac_core(iStart, nCount)

      ! !*-- update core%rho
      ! allocate(self%core%rho(m%cells%numcell_clone))
      ! do iCell = 1, m%cells%numcell_clone
      !    self%core%rho(iCell) = sum(self%frac_core%mass(iCell,:)) / m%cells%vcell(iCell)
      ! end do
      
      ! Allocate and randomize velocities
      allocate(self%core%cell_velocity(m%cells%numcell_clone, m%sim%numvel))
      call random_number(self%core%cell_velocity)
      
      ! Allocate and set deriv_velocity to zero
      allocate(self%core%deriv_velocity(m%cells%numcell_clone, m%sim%numdim, m%sim%numvel))
      self%core%deriv_velocity = 0.0_REAL64
      
      call clone_barrier()
      if (myID == 0) then
         write(*,*) 'Done initializing mesh from file'
      end if
    END ASSOCIATE
  end subroutine init

end module fakemesh

