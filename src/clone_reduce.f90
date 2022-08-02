submodule (clone_lib_module) clone_reduce_module

  use iso_fortran_env, only: INT64, REAL64
#ifdef ENABLE_MPI
  use mpi
#endif
  implicit none

  
contains

  module procedure clone_reduce_i_0
    implicit none
    integer :: my_mpi_Op
    logical :: bcast
    integer :: ierror

    if (present(do_bcast)) then
       bcast = do_bcast
    else
       bcast = .false.
    end if
    
#ifndef ENABLE_MPI
    result_out = value_in
#else
    select case (op)
    case (CLONE_SUM)
       my_mpi_op = MPI_SUM
    case (CLONE_MAX)
       my_mpi_op = MPI_MAX
    case (CLONE_MIN)
          my_mpi_op = MPI_MIN
       case default
          write(*,*) 'Invalid OP sent to clone_reduce'
          return
    end select
    
    if (bcast) then
       call mpi_Allreduce(value_in, result_out, 1, MPI_INTEGER, my_mpi_Op, myComm, ierror)
    else
       call mpi_reduce(value_in, result_out, 1, MPI_INTEGER, my_mpi_Op, 0, myComm, ierror)
    end if
#endif
  end procedure clone_reduce_i_0

  module procedure clone_reduce_i_1
    implicit none
    integer :: tmp_val
    
    select case (op)
    case (CLONE_SUM)
       tmp_val = sum(value_in)
    case (CLONE_MAX)
       tmp_val = maxval(value_in)
    case (CLONE_MIN)
       tmp_val = minval(value_in)
       case default
          write(*,*) 'Invalid OP sent to clone_reduce'
          return
    end select

    call clone_reduce_i_0(result_out, tmp_val, op, do_bcast)
  end procedure clone_reduce_i_1

  module procedure clone_reduce_i64_0
    implicit none
    integer :: my_mpi_Op
    logical :: bcast
    integer :: ierror

#ifndef ENABLE_MPI
    result_out = value_in
#else
    if (present(do_bcast)) then
       bcast = do_bcast
    else
       bcast = .false.
    end if
    
    select case (op)
    case (CLONE_SUM)
       my_mpi_op = MPI_SUM
    case (CLONE_MAX)
       my_mpi_op = MPI_MAX
    case (CLONE_MIN)
          my_mpi_op = MPI_MIN
       case default
          write(*,*) 'Invalid OP sent to clone_reduce'
          return
    end select
    
    if (bcast) then
       call mpi_Allreduce(value_in, result_out, 1, MPI_INTEGER8, my_mpi_Op, myComm, ierror)
    else
       call mpi_reduce(value_in, result_out, 1, MPI_INTEGER8, my_mpi_Op, 0, myComm, ierror)
    end if
#endif
  end procedure clone_reduce_i64_0

  module procedure clone_reduce_i64_1
    implicit none
    integer(INT64) :: tmp_val
    
    select case (op)
    case (CLONE_SUM)
       tmp_val = sum(value_in)
    case (CLONE_MAX)
       tmp_val = maxval(value_in)
    case (CLONE_MIN)
       tmp_val = minval(value_in)
       case default
          write(*,*) 'Invalid OP sent to clone_reduce'
          return
    end select

    call clone_reduce_i64_0(result_out, tmp_val, op, do_bcast)
  end procedure clone_reduce_i64_1
  
  module procedure clone_reduce_i_i64_0
    implicit none
    integer :: my_mpi_Op
    logical :: bcast
    integer :: ierror
    integer(INT64) :: tmp_val
    tmp_val = value_in
    call clone_reduce_i64_0(result_out, tmp_val, op, do_bcast)
  end procedure clone_reduce_i_i64_0

  module procedure clone_reduce_i_i64_1
    implicit none
    integer(INT64) :: tmp_val
    
    select case (op)
    case (CLONE_SUM)
       tmp_val = sum(value_in)
    case (CLONE_MAX)
       tmp_val = maxval(value_in)
    case (CLONE_MIN)
       tmp_val = minval(value_in)
       case default
          write(*,*) 'Invalid OP sent to clone_reduce'
          return
    end select

    call clone_reduce_i64_0(result_out, tmp_val, op, do_bcast)
  end procedure clone_reduce_i_i64_1


  module procedure clone_reduce_r64_0
    implicit none
    integer :: my_mpi_Op
    logical :: bcast
    integer :: ierror
    real(REAL64) :: local_result

#ifndef ENABLE_MPI
    result_out = value_in
#else
    if (present(do_bcast)) then
       bcast = do_bcast
    else
       bcast = .false.
    end if
    
    select case (op)
    case (CLONE_SUM)
       my_mpi_op = MPI_SUM
    case (CLONE_MAX)
       my_mpi_op = MPI_MAX
    case (CLONE_MIN)
          my_mpi_op = MPI_MIN
       case default
          write(*,*) 'Invalid OP sent to clone_reduce'
          return
    end select
    
    if (bcast) then
       call mpi_Allreduce(value_in, result_out, 1, MPI_REAL8, my_mpi_Op, myComm, ierror)
    else
       call mpi_reduce(value_in, result_out, 1, MPI_REAL8, MPI_SUM, 0, myComm, ierror)
    end if
#endif
  end procedure clone_reduce_r64_0

  module procedure clone_reduce_r64_1
    implicit none
    real(REAL64) :: tmp_val
    
    select case (op)
    case (CLONE_SUM)
       tmp_val = sum(value_in)
    case (CLONE_MAX)
       tmp_val = maxval(value_in)
    case (CLONE_MIN)
       tmp_val = minval(value_in)
       case default
          write(*,*) 'Invalid OP sent to clone_reduce'
          return
    end select

    call clone_reduce_r64_0(result_out, tmp_val, op, do_bcast)
  end procedure clone_reduce_r64_1
  

end submodule clone_reduce_module
