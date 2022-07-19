! This is where MPI communications go, for now only serial
module clone_lib_module
  use define_kind
  use var_wrapper_class, only : var_wrapper
  use iso_fortran_env, only: INT64, REAL64
#ifdef EP_MPI
  ! include 'mpif.h'
  use mpi
  use iso_c_binding, only: c_loc
#endif
  implicit none

  

  
  private
  public clone_myid
  public clone_nprocs
  public clone_get
  public clone_init
  public clone_base_init
  public clone_exit
  public clone_barrier

  integer, parameter :: IDLE = 1
  integer, parameter :: SENT = 2
  integer, parameter :: WAITING = 4

#ifdef EP_MPI
  integer :: myComm = MPI_COMM_WORLD
#endif

  integer :: g_nprocs = 1
  integer :: g_myid = 0
  
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
    
#ifndef EP_MPI
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
#ifdef EP_MPI
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
#ifdef EP_MPI
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
    integer(c_int64_t) :: id_nbr, iClone
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

    iClone = 0
    do iDim=1,m%sim%numdim
       do iTmp=1,m%levels%numtop
          
          iCell = m%levels%ltop(iTmp)
          ! Low side
          id_nbr = nbrs(iCell, 2 * iDim - 1)
          if (id_nbr < iStart .or. id_nbr > iEnd) then
             ! Off processor
             call update_clone_id(id_nbr, iClone, clone_map)
             nbrs(iCell, 2 * iDim - 1) = id_nbr + m%cells%numcell
          else
             ! On processor
             nbrs(iCell, 2 * iDim - 1) = id_nbr - iStart + 1
          end if

          ! High side
          id_nbr = nbrs(iCell, 2 * iDim)
          if (id_nbr < iStart .or. id_nbr > iEnd) then
             ! Off processor
             call update_clone_id(id_nbr, iClone, clone_map)
             nbrs(iCell, 2 * iDim ) = id_nbr + m%cells%numcell
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
    subroutine update_clone_id(the_ID, the_count, the_clone_map)
      ! Checks through the clone map to see if the_ID exists already.
      ! If not, adds it to the_clone_map, increments the_count
      use iso_fortran_env, only: INT64
      implicit none

      integer(INT64), intent(inout) :: the_ID
      integer(INT64), intent(inout) :: the_count
      integer(INT64), intent(inout) :: the_clone_map(:)
      
      integer :: j
      integer(INT64) :: new_clone_id

      new_clone_id = -1

      ! Search backwards through the clone array
      ! Must be a smarter way to do this.
      ! This is an n^2 search.  
      do j = the_count, 1, -1
         if (the_clone_map(j) == the_ID) then
            new_clone_id = j
            exit
         end if
      end do

      ! Increment the_count if we need to
      if (new_clone_id < 0) then
         the_count = the_count + 1
         the_clone_map(the_count) = the_ID
         new_clone_id = the_count
      end if

      the_ID = new_clone_id

    end subroutine update_clone_id
    integer(INT64) function get_clone_id(the_ID, the_count, the_clone_map)
      ! Checks through the clone map to see if the_ID exists already.
      ! If not, adds it to the_clone_map, increments the_count
      use iso_fortran_env, only: INT64
      implicit none

      integer(INT64), intent(in) :: the_ID
      integer(INT64), intent(inout) :: the_count
      integer(INT64), intent(inout) :: the_clone_map(:)

      integer :: j

      get_clone_id = -1

      ! Search backwards through the clone array
      ! Must be a smarter way to do this.
      ! This is an n^2 search.  
      do j = the_count, 1, -1
         if (the_clone_map(j) == the_ID) then
            get_clone_id = j
            exit
         end if
      end do

      ! Increment the_count if we need to
      if (get_clone_id < 0) then
         the_count = the_count + 1
         the_clone_map(the_count) = the_ID
         get_clone_id = the_count
      end if

    end function get_clone_id

  end subroutine update_nbrs_and_get_clone_map
  
  subroutine clone_init(m, nbrs, myid, partition)
    ! Initializes the clone arrays and
    ! fixes the neighboring cell IDs 
    use iso_fortran_env, only: INT64
    use mesh_types, only: mesh_t
    implicit none
    type(mesh_t), intent(inout) :: m
    integer(INT64), allocatable :: nbrs(:,:)
    integer(INT64) :: partition(0:)
    integer, intent(in) :: myid
    
    
    integer :: l, nprocs, ierror
    integer :: iDim, iProc, iNode, iTmp, iNow, index
    integer(INT64) :: id_lo, id_hi, iCell, iStart, iEnd
    
    integer(INT64), allocatable :: clone_map(:)
    integer, allocatable :: proc_map(:), tmp_recv(:), tmp_id_recv(:)
    
    ASSOCIATE(                                   &
         numtop => m%levels%numtop,              &
         numcell_clone => m%cells%numcell_clone, &
         numcell => m%cells%numcell              &
         )

      ! Initialize convenience scalars
      nprocs = size(partition, 1)
      write(*,*) 'myid=', myid, ': partition=', partition
      iStart = partition(myid)
      iEnd = partition(myid + 1) - 1

      ! Update the neighbor array and get a clone map of
      ! clone cell ID (from m%cells%numcell +1 -> m%cells%numcellclone)
      ! to absolute cell number
      call update_nbrs_and_get_clone_map(m, iStart, iEnd, nbrs, clone_map)

      ! Count clones by proc
      allocate(tmp_recv(0:nprocs))
      n_nodes = 0
      tmp_recv = 0
      do iDim=1, m%sim%numdim
         do iTmp = numcell + 1, numcell_clone
            iCell = clone_map(iTmp - numcell)
            iProc = get_proc_id(iCell, nprocs, partition)
            if ( iProc /= myid) then
               if (tmp_recv(iProc) == 0) then
                  ! only check tmp_recv because recv and send procs
                  ! are same due to geometric constraint
                  n_nodes = n_nodes + 1
               end if
               tmp_recv(iProc) = tmp_recv(iProc) + 1
            end if
         end do
      end do

      ! Generate the Node structure
      ! Replaces tmp with a mapping to node
      ! Repurposes tmp_recv as index
      allocate(nodes(n_nodes))
      allocate(proc_map(0:nprocs-1))
      proc_map = -1
      iNow = 0
      iNode = 0
      do iDim=1, m%sim%numdim
         do iTmp = numcell + 1, numcell_clone
            iCell = clone_map(iTmp-numcell)
            iProc = get_proc_id(iCell, nprocs, partition)
            iNode = proc_map(iProc)
            if (iNode < 0) then
               ! initialize node structure and repurpose tmp_recv(iProc)
               iNow = iNow + 1
               iNode = iNow
               proc_map(iProc) = iNode
               nodes(iNode)%rank = iProc
               nodes(iNode)%status = IDLE
               nodes(iNode)%nRecv = tmp_recv(iProc)
               allocate(tmp_id_recv(tmp_recv(iProc)))
               allocate(nodes(iNode)%recv_map(tmp_recv(iProc)))
               tmp_recv(iProc) = 0
            end if
            tmp_recv(iProc) = tmp_recv(iProc) + 1
            index = tmp_recv(iProc)
            tmp_id_recv(index) = iCell - partition(iProc) + 1
            nodes(iNode)%recv_map(index) = iTmp
         end do
      end do

      ! Deallocate temporary memory
      deallocate(clone_map)
      deallocate(proc_map)

#ifdef EP_MPI
      
      ! Now Ask other processors what to send
      ! and let them know what we expect to receive
      ! repurpose tmp_recv as a MPI request
      
      do iNode = 1, n_nodes
         iProc = nodes(iNode)%rank
         nodes(iNode)%nRecv = tmp_recv(iProc)
         ! post receive for how many we need to send
         call mpi_irecv(nodes(iNode)%nSend, 1, MPI_INTEGER, &
             nodes(iNode)%rank, 1, myComm, nodes(iNode)%request_recv, ierror)
         ! post send for how many we expect to receive
         call mpi_isend(nodes(iNode)%nRecv, 1, MPI_INTEGER, &
              nodes(iNode)%rank, 1, myComm, tmp_recv(iNode), ierror)
         ! post IDs of cells we expect to receive
         call mpi_isend(tmp_id_recv, nodes(iNode)%nRecv, MPI_INTEGER, &
             nodes(iNode)%rank, 3, myComm, nodes(iNode)%request_send, ierror)
      end do

      call flush()
      call clone_barrier()
      ! Wait for communications to end and collect up what we need to send
      do iNode=1, n_nodes
         iProc = nodes(iNode)%rank
         call mpi_wait(nodes(iNode)%request_recv, MPI_STATUS_IGNORE, ierror)
         
         ! Allocate space and post receive for IDS
         allocate(nodes(iNode)%send_id(nodes(iNode)%nSend))
         call mpi_irecv(nodes(iNode)%send_id, nodes(iNode)%nSend, MPI_INTEGER, &
              nodes(iNode)%rank, 3, myComm, nodes(iNode)%request_recv, ierror)
      end do

      ! Wait for all communications to finish
      call mpi_waitall(n_nodes, tmp_recv(1:n_nodes), MPI_STATUSES_IGNORE, ierror)
      call mpi_waitall(n_nodes, nodes(:)%request_recv, MPI_STATUSES_IGNORE, ierror)
      call mpi_waitall(n_nodes, nodes(:)%request_send, MPI_STATUSES_IGNORE, ierror)
      call flush()
      call clone_barrier()
#endif

      ! Deallocate temporary data structures
      if ( allocated(tmp_id_recv) ) deallocate(tmp_id_recv)
      if ( allocated(tmp_recv) ) deallocate(tmp_recv)
         
    END ASSOCIATE

  CONTAINS
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
#ifdef EP_MPI
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
#ifdef EP_MPI
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
#ifdef EP_MPI
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

