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

! This is where MPI communications go, for now only serial
module clone_lib_module
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
  public clone_abort

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

  subroutine clone_abort(message)
    implicit none
    integer :: ierr
    character(len=*) :: message
    write(*,*) g_myid, '_____________ERROR ERROR ERROR________________'
    write(*,*) g_myid, message
    write(*,*) g_myid, '----------------------------------------------'
    call flush(6)
#ifdef ENABLE_MPI
    call MPI_Abort(mycomm, 1, ierr)
#else
    stop 'Error!  See above.'
#endif
  end subroutine clone_abort
  
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
       stop 'wrong data type'
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
    if (allocated(nodes)) then
       do i = 1, n_nodes
          if (allocated(nodes(i)%send_id)) deallocate(nodes(i)%send_id)
          if (allocated(nodes(i)%recv_map)) deallocate(nodes(i)%recv_map)
       end do
       deallocate(nodes)
    end if
    return
  end subroutine clone_exit

  subroutine clone_barrier()
    implicit none
#ifdef ENABLE_MPI
    integer :: ierror
    call MPI_Barrier(myComm, ierror)
#endif
  end subroutine clone_barrier
    
  subroutine clone_init(m, nbrs, clone_map)
    ! Initializes the clone arrays and
    ! fixes the neighboring cell IDs 
    use iso_fortran_env, only: INT64
    use mesh_types, only: mesh_t
    use binreader, only: binfile
    implicit none
    type(mesh_t), intent(inout) :: m
    integer(INT64), pointer :: nbrs(:,:)
    integer(INT64), pointer :: clone_map(:)
    type(binFile) :: bfp
    integer :: l, ierror
    integer :: iDim, iProc, iNode, iTmp, iNow, index
    integer(INT64) :: id_lo, id_hi, iCell, iStart, iEnd
    
    integer, allocatable :: proc_map(:), tmp(:)
    type(data_t), allocatable :: tmp_id_recv(:)
    integer, parameter :: TAG_NSEND=1, TAG_IDS=3
    
    ASSOCIATE(                                   &
         numtop => m%levels%numtop,              &
         allnumtop => m%levels%allnumtop,        &
         numcell_clone => m%cells%numcell_clone, &
         numcell => m%cells%numcell,             &
         partition => m%cells%cell_address,      &
         nprocs => g_nprocs                      &
         )

      call clone_barrier()

      ! Initialize convenience scalars
      iStart = partition(g_myid)
      iEnd = partition(g_myid + 1) - 1

      ! Next two loops loop over remote cells, so we don't have to
      ! check if they are on-processor

      ! Count clones by proc
      allocate(tmp(0:nprocs-1))
      n_nodes = 0
      tmp = 0
      do iTmp = 1, allnumtop - numtop
         iCell = clone_map(iTmp)
         iProc = get_proc_id(iCell, nprocs, partition)
         if (iProc < 0) then
            write(*,*) 'idx=',iTmp
            write(*,*) 'map=', clone_map
            write(*,*)  'negative iProc!', iCell, iProc
            call flush(6)
            call clone_abort('negative iProc')
         end if
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
      ! write(*,*) g_myid, 'clonemap=', clone_map
      ! write(*,*) g_myid, 'sizeof_clonemap', shape(clone_map)
      ! write(*,*) g_myid, 'numbers:', allnumtop - numtop, numcell_clone - numcell

      do iTmp = numcell + 1, numcell_clone
         iCell = clone_map(iTmp-numcell)               ! Global cell ID of clone
         iProc = get_proc_id(iCell, nprocs, partition) ! Remote processor ID of clone
         if (iProc == g_myID) then
            write(*,'("Problem? ONPE clone: ", i6, 3(i8,","))') iProc, iCell, iStart, iEnd
            cycle
         end if
         if (proc_map(iProc) < 0) then
            ! First neighbor for given processor:
            ! initialize node structure and repurpose tmp(iProc)
            iNow = iNow + 1
            proc_map(iProc) = iNow
            nodes(iNow)%rank = iProc
            nodes(iNow)%status = IDLE
            nodes(iNow)%nRecv = tmp(iProc)
            
            ! Allocate space for mapping data received
            if (allocated(nodes(iNow)%recv_map)) then
               deallocate(nodes(iNow)%recv_map)
            end if
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
         call clone_abort('unequal nodes')
      end if
      
      ! Deallocate temporary memory
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

