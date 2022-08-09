! This is where MPI communications go, for now only serial
module clone_lib_module
  
  use define_kind
  use var_wrapper_class, only : var_wrapper
  use iso_fortran_env, only: INT64, REAL64
#ifdef ENABLE_MPI
  ! include 'mpif.h'
  use mpi
  use iso_c_binding, only: c_loc
#endif
  
  
  public clone_myid
  public clone_nprocs
  public clone_get
  public clone_init
  public clone_base_init
  public clone_exit
  public clone_barrier
  public clone_reduce

  integer, parameter :: IDLE = 1
  integer, parameter :: SENT = 2
  integer, parameter :: WAITING = 4


#ifdef ENABLE_MPI
  integer :: myComm = MPI_COMM_WORLD
#else
  integer :: myComm = -1
#endif

  integer :: g_nprocs = 1
  integer :: g_myid = 0

  ! Reductions
  integer, parameter :: CLONE_SUM = 1
  integer, parameter :: CLONE_MAX = 2
  integer, parameter :: CLONE_MIN = 3
  interface clone_reduce
     module procedure clone_reduce_i_0
     module procedure clone_reduce_i_1
     module procedure clone_reduce_i_i64_0
     module procedure clone_reduce_i_i64_1
     module procedure clone_reduce_i64_0
     module procedure clone_reduce_i64_1
     module procedure clone_reduce_r64_0
     module procedure clone_reduce_r64_1
  end interface clone_reduce
  
  
  ! Gets
  interface clone_get
     procedure clone_get_vw
     procedure clone_get_2
     
     procedure clone_get_i32_1
     procedure clone_get_i64_1
     
     procedure clone_get_r64_1
  end interface clone_get

  ! Data structure to hold remote data
  integer, parameter :: INDEX_SEND = 1
  integer, parameter :: INDEX_RECV = 2
  
  integer, parameter :: DATA_I = 1
  integer, parameter :: DATA_I64 = 2
  integer, parameter :: DATA_R64 = 4
  type :: data_t
     ! Data for transferring over MPI
     ! I could create an abstract type
     ! with specializations for the different
     ! types, but that would be overkill for
     ! something this simple
     integer :: n
     integer, pointer, dimension(:) :: i => NULL()
     integer(INT64), pointer, dimension(:) :: i64 => NULL()
     real(REAL64), pointer, dimension(:) :: r64 => NULL()
     contains
       procedure :: alloc => data_alloc
       procedure :: release => data_release
  end type data_t
  
  type :: node_t
     integer :: rank                       ! Remote rank
     integer :: status                   ! Unused, must remove
     integer :: nSend                    ! How many to send
     integer :: nRecv                    ! How many to receive
     integer :: request_send             ! send mpi request ID
     integer :: request_recv             ! recv mpi request ID
     integer, allocatable :: send_id(:)  ! Ids of cell sto send
     integer, allocatable :: recv_map(:) ! Where received cells fit into my data
  end type node_t

  integer :: n_nodes                             ! Number of neighboring PEs
  type(node_t), target, allocatable, dimension(:) :: nodes ! space for sending / receiving data


  interface
     module subroutine clone_reduce_i_0(result_out, value_in, op, do_bcast)
       implicit none
       integer, intent(out) :: result_out
       integer, intent(in) :: value_in
       integer, intent(in) :: op
       logical, optional, intent(in) :: do_bcast
     end subroutine clone_reduce_i_0
     module subroutine clone_reduce_i_1(result_out, value_in, op, do_bcast)
       implicit none
       integer, intent(out) :: result_out
       integer, intent(in) :: value_in(:)
       integer, intent(in) :: op
       logical, optional, intent(in) :: do_bcast
     end subroutine clone_reduce_i_1

     module subroutine clone_reduce_i_i64_0(result_out, value_in, op, do_bcast)
       implicit none
       integer(INT64), intent(out) :: result_out
       integer, intent(in) :: value_in
       integer, intent(in) :: op
       logical, optional, intent(in) :: do_bcast
     end subroutine clone_reduce_i_i64_0
     module subroutine clone_reduce_i_i64_1(result_out, value_in, op, do_bcast)
       implicit none
       integer(INT64), intent(out) :: result_out
       integer, intent(in) :: value_in(:)
       integer, intent(in) :: op
       logical, optional, intent(in) :: do_bcast
     end subroutine clone_reduce_i_i64_1

     module subroutine clone_reduce_i64_0(result_out, value_in, op, do_bcast)
       use iso_fortran_env, only: INT64
       implicit none
       integer(INT64), intent(out) :: result_out
       integer(INT64), intent(in) :: value_in
       integer, intent(in) :: op
       logical, optional, intent(in) :: do_bcast
     end subroutine clone_reduce_i64_0
     module subroutine clone_reduce_i64_1(result_out, value_in, op, do_bcast)
       implicit none
       integer(INT64), intent(out) :: result_out
       integer(INT64), intent(in) :: value_in(:)
       integer, intent(in) :: op
       logical, optional, intent(in) :: do_bcast
     end subroutine clone_reduce_i64_1

     module subroutine clone_reduce_r64_0(result_out, value_in, op, do_bcast)
       use iso_fortran_env, only: REAL64
       implicit none
       real(REAL64), intent(out) :: result_out
       real(REAL64), intent(in) :: value_in
       integer, intent(in) :: op
       logical, optional, intent(in) :: do_bcast
     end subroutine clone_reduce_r64_0
     module subroutine clone_reduce_r64_1(result_out, value_in, op, do_bcast)
       implicit none
       real(REAL64), intent(out) :: result_out
       real(REAL64), intent(in) :: value_in(:)
       integer, intent(in) :: op
       logical, optional, intent(in) :: do_bcast
     end subroutine clone_reduce_r64_1
  end interface

contains
  
  integer function clone_nprocs()
    implicit none
    clone_nprocs = g_nprocs
  end function clone_nprocs

  integer function clone_myid()
    implicit none
    clone_myid = g_myid
  end function clone_myid

  function data_array_alloc(my_type) result(ptr)
    ! allocates n_nodes arrays of right size
    ! for MPI transfers.  This should really
    ! not be necessary, but the Fortran MPI
    ! interface does not allow for array slices
    ! to be passed in
    implicit none
    type(data_t), pointer, dimension(:,:) :: ptr
    integer, intent(in) :: my_type

    integer :: i
    allocate(ptr(n_nodes, 2))
    do i = 1, n_nodes
       call ptr(i,INDEX_SEND)%alloc(nodes(i)%nSend, my_type)
       call ptr(i,INDEX_RECV)%alloc(nodes(i)%nRecv, my_type)
    end do
  end function data_array_alloc
    
  subroutine data_array_dealloc(ptr)
    ! deallocates the array of data
    implicit none
    type(data_t), pointer, dimension(:,:) :: ptr
    integer :: isz
    integer :: i
    isz = size(ptr, 1)
    do i = 1, isz
       call ptr(i,INDEX_SEND)%release()
       call ptr(i,INDEX_RECV)%release()
    end do
    deallocate(ptr)
  end subroutine data_array_dealloc
    
  subroutine data_alloc(self, n, the_type)
    ! allocates a single data element
    ! of given size and type

    implicit none
    class(data_t) :: self
    integer, intent(in) :: n, the_type
    self%n = n

    nullify(self%i)
    nullify(self%i64)
    nullify(self%r64)
    if (the_type == DATA_I) then
       allocate(self%i(n))
    else if (the_type == DATA_I64) then
       allocate(self%i64(n))
    else if (the_type == DATA_R64) then
       allocate(self%r64(n))
    else
       write(*,*) 'wrong data type'
    end if
  end subroutine data_alloc
  
  subroutine data_release(self)
    class(data_t) :: self
    if (associated(self%i)) deallocate(self%i)
    if (associated(self%i64)) deallocate(self%i64)
    if (associated(self%r64)) deallocate(self%r64)
    nullify(self%i)
    nullify(self%i64)
    nullify(self%r64)
  end subroutine data_release

  pure integer function get_proc_id(theCell, nprocs, partition)
    use iso_fortran_env, only: INT64
    implicit none
    integer(INT64), intent(in) :: theCell
    integer, intent(in) :: nprocs
    integer(INT64), intent(in) :: partition(0:)
    integer :: i

    get_proc_id = -1
    do i = 1, nprocs
       if ( theCell < partition(i) ) then
          get_proc_id = i - 1
          exit
       end if
    end do
  end function get_proc_id
    
#ifndef ENABLE_MPI
  subroutine clone_base_init(myid, nprocs)
    ! no MPI, so no work
    implicit none
    integer, intent(inout) :: myid, nprocs
    myid = 0
    nprocs = 1
    return
  end subroutine clone_base_init
#else
  subroutine clone_base_init(myid, nprocs)
    ! Initializes MPI if needed and returns myid and nprocs
    implicit none
    integer, intent(out) :: myid, nprocs
    integer :: ierror = 0
    call MPI_INIT(ierror)
    myComm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(myComm, nprocs, ierror)
    call MPI_COMM_RANK(myComm, myid, ierror)
    g_myid = myid
    g_nprocs = nprocs
  end subroutine clone_base_init
#endif
  subroutine clone_exit()
    implicit none
    integer :: i
#ifdef ENABLE_MPI
    ! close the MPI communication
    call MPI_FINALIZE(i)
#endif
    ! Deallocate data structure
    do i = 1, n_nodes
       if (allocated(nodes(i)%send_id)) deallocate(nodes(i)%send_id)
       if (allocated(nodes(i)%recv_map)) deallocate(nodes(i)%recv_map)
    end do
    deallocate(nodes)
    return
  end subroutine clone_exit

  subroutine clone_barrier()
    implicit none
#ifdef ENABLE_MPI
    integer :: ierror
    call MPI_Barrier(myComm, ierror)
#endif
  end subroutine clone_barrier
    
  subroutine update_nbrs_and_get_clone_map(m, iStart, iEnd, nbrs, clone_map)
    ! Counts m%cells%numcell_clone,
    ! Changes nbrs to be range 1 - m%cells%numcell_clone
    ! Updates m%allnumtop
    ! Allocates and initializes m%alltop
    ! initializes clone_map to map icell to real cell id

    use iso_fortran_env, only: INT64
    use iso_c_binding, only: c_int64_t
    use mesh_types, only: mesh_t
    implicit none
    type(mesh_t), intent(inout) :: m
    integer(c_int64_t), intent(in) :: iStart, iEnd
    integer(INT64), intent(inout) :: nbrs(:,:)
    integer(INT64), allocatable, intent(inout) :: clone_map(:)
    integer(c_int64_t) :: id_nbr, iClone, id_new
    integer :: iCell, nClones, iDim, iTmp
    logical :: found

    ! Get an upper bound on number of clones
    nClones = 0
    do iDim=1,m%sim%numdim
       do iTmp=1,m%levels%numtop
          iCell = m%levels%ltop(iTmp)
          
          ! Low side
          id_nbr = nbrs(iCell, 2 * iDim - 1)
          if (id_nbr < iStart .or. id_nbr > iEnd) then
             nClones = nClones + 1
          end if

          ! High side
          id_nbr = nbrs(iCell, 2 * iDim)
          if (id_nbr < iStart .or. id_nbr > iEnd) then
             nClones = nClones + 1
          end if
       end do
    end do

    ! nClones is now an upper bound on number of clones
    if (allocated(clone_map)) deallocate(clone_map)
    allocate(clone_map(nClones))

    ! This loop converts nbrs array from absolute cell
    ! number to local IDS running from 1 -> numcell_clone
    !
    ! The clone_map array will map local clone IDs (numcell+1 -> numcell_clone)
    ! to global cell IDs primarily for processor identification purposes
    !

    iClone = 0
    do iDim=1,m%sim%numdim
       do iTmp=1,m%levels%numtop
          
          iCell = m%levels%ltop(iTmp)
          
          ! Low side
          id_nbr = nbrs(iCell, 2 * iDim - 1)
          if (id_nbr < iStart .or. id_nbr > iEnd) then
             ! Off processor
             call generate_clone_id(id_nbr, id_new, iClone, clone_map)
             nbrs(iCell, 2 * iDim - 1) = id_new + m%cells%numcell
          else
             ! On processor
             nbrs(iCell, 2 * iDim - 1) = id_nbr - iStart + 1
          end if

          ! High side
          id_nbr = nbrs(iCell, 2 * iDim)
          if (id_nbr < iStart .or. id_nbr > iEnd) then
             ! Off processor
             call generate_clone_id(id_nbr, id_new, iClone, clone_map)
             nbrs(iCell, 2 * iDim ) = id_new + m%cells%numcell
          else
             ! On processor
             nbrs(iCell, 2 * iDim) = id_nbr - iStart + 1
          end if
       end do
    end do

    ! Update numcell_clone, allnumtop, and alltop
    m%cells%numcell_clone = iClone + m%cells%numcell
    m%levels%allnumtop = m%levels%numtop + iClone
    allocate(m%levels%alltop(m%levels%allnumtop))
    m%levels%alltop(1:m%levels%numtop) = m%levels%ltop(1:m%levels%numtop)
    do iCell = 1, iClone
       m%levels%alltop(iCell + m%levels%numtop) = iCell + m%cells%numcell
    end do

  contains
    pure subroutine generate_clone_id(old_clone_id, new_clone_id, the_count, the_clone_map)
      ! Checks through the clone map to see if old_clone_id exists already.
      ! If not, adds it to the_clone_map, increments the_count
      use iso_fortran_env, only: INT64
      implicit none

      integer(INT64), intent(in)    :: old_clone_id
      integer(INT64), intent(out)   :: new_clone_id
      integer(INT64), intent(inout) :: the_count
      integer(INT64), intent(inout) :: the_clone_map(:)
      
      integer :: j

      new_clone_id = -1

      ! Search backwards through the clone array
      ! Must be a smarter way to do this.
      ! This is an n^2 search.  
      do j = the_count, 1, -1
         if (the_clone_map(j) == old_clone_id) then
            new_clone_id = j
            exit
         end if
      end do

      ! Increment the_count if we need to
      if (new_clone_id < 0) then
         the_count = the_count + 1
         the_clone_map(the_count) = old_clone_id
         new_clone_id = the_count
      end if

    end subroutine generate_clone_id
    
  end subroutine update_nbrs_and_get_clone_map

  subroutine clone_update_AMR_boundary(m, nbrs, clone_map, pioid)
    ! Addes additional clones to cells that lie
    ! on the boundary
    use iso_fortran_env, only: INT64
    use mesh_types, only: mesh_t
    use pio_interface, only: pio_get_range_i64, pio_release
    implicit none

    type(mesh_t) :: m
    integer(INT64) :: nbrs(:,:)
    integer(INT64), intent(in) :: clone_map(:)
    integer, intent(in) :: pioid
    
    ! magic offsets for xRage mesh
    integer, parameter :: offsets(3,3) =  reshape([2,4,6, 1,4,5, 1,2,3],[3,3])

    integer :: iTop, iDim, iProc, iBase, iNode, jTmp, iClone, ierror, iLast, iTmp
    integer :: n_shift, n_additional, n_external, old_id
    integer(INT64) :: iCell, iNbr
    integer(INT64), pointer :: tmp_cell_level(:)
    integer, pointer :: old_cell_level(:)

    integer, allocatable :: clone_to_proc(:), new_index(:)
    integer, allocatable :: node_map(:), node_count(:), sister_clones(:,:)
    type(data_t), allocatable :: new_remote_id(:)
    integer, allocatable :: request_send(:), new_recv_map(:), new_send_size(:)
    integer(INT64), pointer :: daughter(:) => null()

    ! tags for messages
    integer, parameter :: TAG_N=101, TAG_IDS=102

    ASSOCIATE( &
         ndim => m%sim%numdim, &
         numcell => m%cells%numcell, &
         numcell_clone => m%cells%numcell_clone, &
         partition => m%cells%cell_address, &
         numtop => m%levels%numtop, &
         ltop => m%levels%ltop, &
         allnumtop => m%levels%allnumtop, &
         alltop => m%levels%alltop, &
         cell_level => m%levels%cell_level &
         )

      !*-- Update AMR cell_level and insert clones for coarse cells at high boundaries

      ! At this point our clone cells only contain one of the
      ! (1/2/4) T-Cell faces in 1/2/3D. We will now add in the
      ! reamaining faces.
      !
      ! For this we need to know the levels of our clone cells
      
      !* The cell levels will be overridden with final numbers at the end

      ! Read in cell level from the PIO file and populate m%levels%cell_level(1:numcell)
      tmp_cell_level => pio_get_range_i64(pioid, "cell_level", 0, partition(g_myid), &
           partition(g_myid+1)-partition(g_myid))
      allocate(m%levels%cell_level(numcell_clone))
      m%levels%cell_level(1:numcell) = tmp_cell_level
      call pio_release(tmp_cell_level)

      ! Fill in the clone cell levels
      call clone_get(m%levels%cell_level)

      if (ndim < 2 .or. numcell_clone == numcell) then
         ! no changes required
         return
      end if
#ifdef ENABLE_MPI

      ! convenience scalar
      n_external = numcell_clone - numcell
      
      ! Set node counts to existing values and map clones to nodes
      allocate(node_count(n_nodes))
      allocate(node_map(n_external)) 
      node_map = -1
      do iNode = 1, n_nodes
         iProc = nodes(iNode)%rank
         node_count(iNode) = nodes(iNode)%nRecv
         node_map(nodes(iNode)%recv_map - numcell) = iNode
      end do

      ! initialize the index map that will provide new indices of old clones
      allocate(new_index(n_external))
      do iCell = 1, n_external
         new_index(iCell) = iCell + numcell
      end do

      n_additional = 0
      n_shift = 2 * ndim - 3 ! 1 extra in 2D; 3 extra in 3D; 1D was eliminated above
      allocate(sister_clones(n_shift, n_external)) ! IDs of new clone cells
      sister_clones = 0
      ! calculate which clones need sister clones
      do iDim = 1, ndim
         do iTop = 1, numtop

            iCell = ltop(iTop)

            ! Lo side
            iNbr = nbrs(iCell, 2 * iDim - 1) 
            if ((iNbr > numcell) .and. (m%levels%cell_level(iNbr) > m%levels%cell_level(iCell))) then
               iNbr = iNbr - numcell
               iProc = node_map(iNbr)
               node_count(iNode) = node_count(iNode) + n_shift ! new count of cells on processor
               n_additional = n_additional + n_shift           ! Number of new clones
               iBase = clone_map(iNbr) - partition(iProc) + 1
               sister_clones(1:n_shift, iNbr) = iBase + offsets(1:n_shift, iDim)
               ! Clones higher than iNbr have to be shifted
               new_index(iNbr + 1 : n_external) = new_index(iNbr + 1 : n_external) + n_shift
            end if

            ! Hi Side
            iNbr = nbrs(iCell, 2 * iDim)
            if ((iNbr > numcell) .and. (m%levels%cell_level(iNbr) > m%levels%cell_level(iCell))) then
               iNbr = iNbr - numcell
               iProc = node_map(iNbr + numcell)
               node_count(iNode) = node_count(iNode) + n_shift ! new count of cells on processor
               n_additional = n_additional + n_shift           ! Number of new clones
               iBase = clone_map(iNbr) - partition(iProc) + 1
               sister_clones(1:n_shift, iNbr) = iBase + offsets(1:n_shift, iDim)
               ! Clones higher than iNbr have to be shifted
               new_index(iNbr + 1 : n_external) = new_index(iNbr + 1 : n_external) + n_shift 
            end if
         end do
      end do

      ! Deallocate node_map, since no longer needed
      deallocate(node_map)

      ! At this point:
      !  - n_additional = number of new clones needed
      !  - sister_clones = remote IDs of sister clones of given neighbor
      !  - new_index = new ID of existing clones
      !  - node_count holds the new size of data we expect to receive

      ! Now we need to
      !  - reset the neighbor IDs in nbrs based on new_index
      !  - reacalculate what remote IDs we expect to receive from each processor
      !  - send those IDs to the remote processor
      !  - receive the new recv_map from remote processorA

      ! allocate tmp space for MPI communications
      allocate(new_remote_id(n_nodes), request_send(n_nodes), new_send_size(n_nodes))
      request_send = MPI_REQUEST_NULL

      do iNode = 1, n_nodes
         iProc = nodes(iNode)%rank
         ! post send and receive for new size
         call mpi_isend(node_count(iNode), 1, MPI_INTEGER, &
              iProc, TAG_N, myComm, nodes(iNode)%request_send, ierror)
         call mpi_irecv(new_send_size(iNode), 1, MPI_INTEGER, &
              iProc, TAG_N, myComm, nodes(iNode)%request_recv, ierror)

         ! Skip work if no change
         if (node_count(iNode) == nodes(iNode)%nRecv) cycle

         ! compute new remote IDs with clones inserted as required
         ! for type 4 and type 5 faces and communicate
         call new_remote_id(iNode)%alloc(node_count(iNode), DATA_I)
         allocate(new_recv_map(node_count(iNode)))
         iLast = 1
         do iTmp = 1,nodes(iNode)%nRecv
            old_id = nodes(iNode)%recv_map(iTmp)
            iClone = old_id - numcell
            new_recv_map(iLast) = new_index(iClone)
            new_remote_id(iNode)%i(iLast) = clone_map(iClone) - partition(iProc) + 1
            iLast = iLast + 1
            if (sister_clones(1, iClone) > 0) then
               do jTmp = 1, n_shift
                  new_remote_id(iNode)%i(iLast) = sister_clones(jTmp, iClone)
                  new_recv_map(iLast) = new_index(iClone) + jTmp
                  iLast = iLast + 1
               end do
            end if
         end do

         ! send the new_remote_ids to remote processor
         call mpi_isend(new_remote_id(iNode)%i, node_count(iNode), MPI_INTEGER, &
              iProc, TAG_IDS, myComm, request_send(iNode), ierror)

         ! Update receive counts and receive map
         nodes(iNode)%nrecv = iLast
         call move_alloc(new_recv_map, nodes(iNode)%recv_map)

      end do

      ! Deallocate sister_clones, since no longer needed
      deallocate(sister_clones)

      ! Now receive the data from the remote and post request for send IDs
      do iNode = 1, n_nodes
         iProc = nodes(iNode)%rank
         call mpi_wait(nodes(iNode)%request_recv, MPI_STATUS_IGNORE, ierror)
         if (new_send_size(iNode) == nodes(iNode)%nSend) then
            nodes(iNode)%request_recv = MPI_REQUEST_NULL
            cycle
         end if
         nodes(iNode)%nsend = new_send_size(iNode)
         deallocate(nodes(iNode)%send_id)
         allocate(nodes(iNode)%send_id(new_send_size(iNode)))
         call mpi_irecv(nodes(iNode)%send_id, nodes(iNode)%nSend, MPI_INTEGER, &
              iProc, TAG_IDS, myComm, nodes(iNode)%request_recv, ierror)
      end do

      ! Update mesh data structures
      numcell_clone = numcell_clone + n_additional
      allnumtop = allnumtop + n_additional
      deallocate(m%levels%alltop)
      allocate(m%levels%alltop(allnumtop))
      m%levels%alltop(1:numtop) = m%levels%ltop(1:numtop)
      do iTmp = 1, (allnumtop-numtop)
         m%levels%alltop(iTmp) = numcell + iTmp
      end do

      ! Now need to update the neighbors array.
      do iTmp = 1, numtop
         iCell = ltop(iTmp)
         do iDim = 1,2*nDim
            iNbr = nbrs(iCell, iDim)
            if (iNbr > numcell) then
               nbrs(iCell,iDim) = new_index(iNbr - numcell)
            end if
         end do
      end do

      ! Deallocate new_index, since no longer needed
      deallocate(new_index)

      ! Wait for MPI receives to finish
      call mpi_waitall(n_nodes, nodes(:)%request_recv, MPI_STATUSES_IGNORE, ierror)

      ! Wait for MPI sends to finish
      call mpi_waitall(n_nodes, nodes(:)%request_send, MPI_STATUSES_IGNORE, ierror)
      call mpi_waitall(n_nodes, request_send, MPI_STATUSES_IGNORE, ierror)


      !*-- Resize AMR cell_level with new ghosts and update clones
      old_cell_level => m%levels%cell_level
      nullify(m%levels%cell_level)
      allocate(m%levels%cell_level(numcell_clone))
      m%levels%cell_level(1:numcell) = old_cell_level(1:numcell)
      deallocate(old_cell_level)
      call clone_get(m%levels%cell_level)

      !*-- Resize daughters to check if any off-processor cells need to be modified
      daughter => m%levels%cell_daughter
      nullify(m%levels%cell_daughter)
      allocate(m%levels%cell_daughter(m%cells%numcell_clone))
      m%levels%cell_daughter(1:m%cells%numcell) = daughter(1:m%cells%numcell)
      call pio_release(daughter)
      call clone_get(m%levels%cell_daughter)

      ! Deallocate memory that was used for MPI buffers
      do iNode = 1, n_nodes
         call new_remote_id(iNode)%release()
      end do
      deallocate(new_remote_id)
      deallocate(node_count)

      ! BLOCK
      !   ! Quick check
      !   use iso_fortran_env, only: INT64
      !   integer(INT64), pointer, dimension(:) :: daughter
      !   daughter => pio_get_range_i64(1, "cell_daughter", 0, &
      !        partition(clone_myid()), partition(clone_myid()+1) - partition(clone_myid()))
      !   do iNode = 1, n_nodes
      !      if (any(daughter(nodes(iNode)%send_id) > 0)) then
      !         write(*,*) 'ERROR ON SEND from ',clone_myid(), nodes(iNode)%rank
      !      end if
      !   end do
      ! END BLOCK

#endif
    END ASSOCIATE

  end subroutine clone_update_AMR_boundary

  subroutine clone_init(m, nbrs, pioid)
    ! Initializes the clone arrays and
    ! fixes the neighboring cell IDs 
    use iso_fortran_env, only: INT64
    use mesh_types, only: mesh_t
    implicit none
    type(mesh_t), intent(inout) :: m
    integer(INT64), allocatable :: nbrs(:,:)
    integer, intent(in) :: pioid
    integer :: l, nprocs, ierror
    integer :: iDim, iProc, iNode, iTmp, iNow, index
    integer(INT64) :: id_lo, id_hi, iCell, iStart, iEnd
    
    integer(INT64), allocatable :: clone_map(:)
    integer, allocatable :: proc_map(:), tmp(:)
    type(data_t), allocatable :: tmp_id_recv(:)
    integer, parameter :: TAG_NSEND=1, TAG_IDS=3
    
    ASSOCIATE(                                   &
         numtop => m%levels%numtop,              &
         numcell_clone => m%cells%numcell_clone, &
         numcell => m%cells%numcell,             &
         partition => m%cells%cell_address       &
         )

      ! Initialize convenience scalars
      nprocs = size(partition, 1)
      iStart = partition(g_myid)
      iEnd = partition(g_myid + 1) - 1

      ! Update the neighbor array and get a clone map of
      ! clone cell ID (from m%cells%numcell +1 -> m%cells%numcellclone)
      ! to absolute cell number
      call update_nbrs_and_get_clone_map(m, iStart, iEnd, nbrs, clone_map)

      ! Next two loops loop over remote cells, so we don't have to
      ! check if they are on-processor

      ! Count clones by proc
      allocate(tmp(0:nprocs-1))
      n_nodes = 0
      tmp = 0
      do iTmp = 1, numcell_clone - numcell
         iCell = clone_map(iTmp)
         iProc = get_proc_id(iCell, nprocs, partition)
         if ( iProc /= g_myid) then
            if (tmp(iProc) == 0) then
               n_nodes = n_nodes + 1
            end if
            tmp(iProc) = tmp(iProc) + 1
         end if
      end do

      ! At this point tmp() holds the number we expect to receive from each processor
      
      ! Generate the Node structure
      ! Replaces tmp with a mapping to node
      ! Restarts tmp counting
      allocate(proc_map(0:nprocs-1), stat=ierror)
      allocate(nodes(n_nodes), stat=ierror)
      allocate(tmp_id_recv(n_nodes), stat=ierror)
      if (ierror /= 0) stop 'error allocating in clone_init()'
      proc_map = -1
      iNow = 0
      do iTmp = numcell + 1, numcell_clone
         iCell = clone_map(iTmp-numcell)               ! Global cell ID of clone
         iProc = get_proc_id(iCell, nprocs, partition) ! Remote processor ID of clone
         if (proc_map(iProc) < 0) then
            ! First neighbor for given processor:
            ! initialize node structure and repurpose tmp(iProc)
            iNow = iNow + 1
            proc_map(iProc) = iNow
            nodes(iNow)%rank = iProc
            nodes(iNow)%status = IDLE
            nodes(iNow)%nRecv = tmp(iProc)
            
            ! Allocate space for mapping data received
            allocate(nodes(iNow)%recv_map(nodes(iNow)%nrecv))

            ! Allocate space for remote IDs of cells we expect to
            ! receive that we will send to the remote processor
            call tmp_id_recv(iNow)%alloc(nodes(iNow)%nrecv, DATA_I)
            tmp_id_recv(iNow)%i = -100

            ! Repurpose tmp(iProc) to hold index in the map above
            ! as we fill it.
            tmp(iProc) = 0
         end if

         ! Insert the local ID of iCell on the remote processor
         ! into the node map and increment the index (tmp)
         iNode = proc_map(iProc)
         tmp(iProc) = tmp(iProc) + 1
         index = tmp(iProc)
         tmp_id_recv(iNode)%i(index) = iCell - partition(iProc) + 1
         nodes(iNode)%recv_map(index) = iTmp
      end do
      if (iNow /= n_nodes) then
         write(*,*) g_myid, '__UNEQUAL NNODES__:    ',iNow, n_nodes, proc_map
      end if
      
      ! Deallocate temporary memory
      deallocate(clone_map)
      deallocate(proc_map)

#ifdef ENABLE_MPI
      
      ! Now Ask other processors what to send
      ! and let them know what we expect to receive
      ! repurpose tmp as a MPI request
      do iNode = 1, n_nodes
         iProc = nodes(iNode)%rank
         ! post receive for how many we need to send
         call mpi_irecv(nodes(iNode)%nSend, 1, MPI_INTEGER, &
             nodes(iNode)%rank, TAG_NSEND, myComm, nodes(iNode)%request_recv, ierror)
         
         ! post how many we expect to receive using tmp array to hold request id
         call mpi_isend(nodes(iNode)%nRecv, 1, MPI_INTEGER, &
              nodes(iNode)%rank, TAG_NSEND, myComm, tmp(iNode), ierror)
         
         ! post IDs of cells we expect to receive
         call mpi_isend(tmp_id_recv(iNode)%i, nodes(iNode)%nRecv, MPI_INTEGER, &
             nodes(iNode)%rank, TAG_IDS, myComm, nodes(iNode)%request_send, ierror)
      end do

      ! Wait for communications to end and collect up what we need to send
      do iNode=1, n_nodes
         iProc = nodes(iNode)%rank
         
         ! Wait for nSend request to finish
         call mpi_wait(nodes(iNode)%request_recv, MPI_STATUS_IGNORE, ierror)
         
         ! Allocate space and post receive for IDS
         !SS Replace mother cells with daughters to handle AMR boundaries
         allocate(nodes(iNode)%send_id(nodes(iNode)%nSend))
         call mpi_irecv(nodes(iNode)%send_id, nodes(iNode)%nSend, MPI_INTEGER, &
              nodes(iNode)%rank, TAG_IDS, myComm, nodes(iNode)%request_recv, ierror)
      end do

      ! Wait for all communications to finish
      ! Sends of number we expect to receive
      call mpi_waitall(n_nodes, tmp(1:n_nodes), MPI_STATUSES_IGNORE, ierror)

      ! Reception of IDs we need to send
      call mpi_waitall(n_nodes, nodes(:)%request_recv, MPI_STATUSES_IGNORE, ierror)

      ! Completion of send of IDs requested from remote
      call mpi_waitall(n_nodes, nodes(:)%request_send, MPI_STATUSES_IGNORE, ierror)
#endif

      ! Deallocate temporary data structures
      if ( allocated(tmp) ) deallocate(tmp)
      do iProc=1,n_nodes
         call tmp_id_recv(iProc)%release()
      end do


      ! Now update clone ids for AMR transfers
      call clone_update_AMR_boundary(m, nbrs, clone_map, pioid)
      
    END ASSOCIATE

  end subroutine clone_init


  
  subroutine clone_get_vw(vw, intFlag)
    type(var_wrapper), intent(inout) :: vw(:)
    integer, intent(in) :: intFlag
  end subroutine clone_get_vw
  subroutine clone_get_2(myValue, nClone, nvars)
    real(REAL64), intent(inout) :: myValue(:,:)
    integer, intent(in) :: nClone
    integer, intent(in) :: nvars
  end subroutine clone_get_2
  
  subroutine clone_get_i32_1(myValue)
    implicit none
    integer, target, intent(inout) :: myValue(:)
#ifdef ENABLE_MPI
    integer :: iNode, iProc, ierror
    type(data_t), pointer :: my_data(:,:)
    
    ! Allocate communication buffers
    my_data => data_array_alloc(DATA_I)
    
    ! post background sends
    do iNode = 1, n_nodes
       ASSOCIATE(myNode => nodes(iNode))
         iProc = myNode%rank
         my_data(iNode, INDEX_SEND)%i = myValue(myNode%send_id)
         call mpi_isend(my_data(iNode,INDEX_SEND)%i, myNode%nSend, MPI_INTEGER, &
              nodes(iNode)%rank, 10, myComm, myNode%request_send, ierror)
         call mpi_irecv(my_data(iNode,INDEX_RECV)%i, myNode%nRecv, MPI_INTEGER, &
              nodes(iNode)%rank, 10, myComm, myNode%request_recv, ierror)
       END ASSOCIATE
    end do

    ! Cycle through nodes waiting for receives
    do iNode = 1, n_nodes
       call mpi_wait(nodes(iNode)%request_recv, MPI_STATUS_IGNORE, ierror)
       myValue(nodes(iNode)%recv_map) = my_data(iNode, INDEX_RECV)%i
    end do

    ! wait for sends to complete
    call mpi_waitall(n_nodes, nodes(:)%request_send, MPI_STATUSES_IGNORE, ierror)

    ! release temp data
    call data_array_dealloc(my_data)
#endif
  end subroutine clone_get_i32_1

  subroutine clone_get_i64_1(myValue)
    use iso_fortran_env, only: INT64
    implicit none
    integer(INT64), target, intent(inout) :: myValue(:)
#ifdef ENABLE_MPI
    integer :: iNode, iProc, ierror
    type(data_t), pointer :: my_data(:,:)

    ! Allocate communication buffers
    my_data => data_array_alloc(DATA_I64)
    
    ! post background sends
    do iNode = 1, n_nodes
       ASSOCIATE(myNode => nodes(iNode))
         my_data(iNode, INDEX_SEND)%i64 = myValue(myNode%send_id)
         call mpi_isend(my_data(iNode, INDEX_SEND)%i64, myNode%nSend, MPI_INTEGER8, &
              myNode%rank, 20, myComm, myNode%request_send, ierror)
         call mpi_irecv(my_data(iNode,INDEX_RECV)%i64, myNode%nRecv, MPI_INTEGER8, &
              nodes(iNode)%rank, 20, myComm, myNode%request_recv, ierror)
       END ASSOCIATE
    end do

    ! Post blocking receives
    do iNode = 1, n_nodes
       call mpi_wait(nodes(iNode)%request_recv, MPI_STATUS_IGNORE, ierror)
       myValue(nodes(iNode)%recv_map) = my_data(iNode, INDEX_RECV)%i64
    end do

    ! wait for sends to complete
    call mpi_waitall(n_nodes, nodes(:)%request_send, MPI_STATUSES_IGNORE, ierror)

    ! release temp data
    call data_array_dealloc(my_data)
        
#endif
  end subroutine clone_get_i64_1
  
  subroutine clone_get_r64_1(myValue)
    use iso_fortran_env, only: REAL64
    implicit none
    real(real64), target, intent(inout) :: myValue(:)
#ifdef ENABLE_MPI
    integer :: iNode, iProc, ierror, nSend, nRecv
    type(data_t), pointer :: my_data(:,:)

    ! Allocate communication buffers
    my_data => data_array_alloc(DATA_R64)
    
    ! post background sends
    do iNode = 1, n_nodes
       ASSOCIATE(myNode => nodes(iNode))
         iProc = myNode%rank
         my_data(iNode,INDEX_SEND)%r64 = myValue(myNode%send_id)
         call mpi_isend(my_data(iNode, INDEX_SEND)%r64, myNode%nSend, MPI_REAL8, &
              myNode%rank, 30, myComm, myNode%request_send, ierror)
         call mpi_irecv(my_data(iNode,INDEX_RECV)%r64, myNode%nRecv, MPI_REAL8, &
              nodes(iNode)%rank, 30, myComm, myNode%request_recv, ierror)
       END ASSOCIATE
    end do

    ! Cycle through nodes waiting for receives
    do iNode = 1, n_nodes
       call mpi_wait(nodes(iNode)%request_recv, MPI_STATUS_IGNORE, ierror)
       myValue(nodes(iNode)%recv_map) = my_data(iNode, INDEX_RECV)%r64
    end do

    ! wait for sends to complete
    call mpi_waitall(n_nodes, nodes(:)%request_send, MPI_STATUSES_IGNORE, ierror)

    ! release temp data
    call data_array_dealloc(my_data)
    
#endif
  end subroutine clone_get_r64_1
end module clone_lib_module

