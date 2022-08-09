######################################
# CMakeLists.txt for EAP Patterns
######################################
cmake_minimum_required(VERSION 3.16.0)
project(mygrad LANGUAGES Fortran CXX C)

option(EP_PIO "Build targets with PIO integration;" "ON")
option(EP_MPI "Build targets with MPI integration;" "ON")
option(EP_SANITIZE "Add address sanitizer flag;" "ON")
option(EP_GZ "Enable GZ file reading in PIO;" "ON")

# unimplemented options
option(EP_HDF5 "Build targets with HDF5 integration; not implemented yet" "OFF")
option(EP_OPENMP "Build targets with OpenMP integration; not implemented yet" "OFF")

# Specify default build
set(CMAKE_BUILD_TYPE Release CACHE STRING "Buld type")

# Specify sources
set(F90_SOURCES
  src/clone_lib.f90
  src/clone_reduce.f90
  src/define_kind.f90 
  src/derivatives_common_template.f90
  src/fakemesh.f90 
  src/main.F90
  src/mem_release.f90
  src/mesh_types.f90
  src/tests.f90
  src/util.f90
  )

set_source_files_properties(
  file1.fpp file2.fpp
  PROPERTIES Fortran_PREPROCESS ON
)
set_source_files_properties(
  ${F90_SOURCES}
  PROPERTIES Fortran_PREPROCESS ON
)

if (EP_SANITIZE )
  if ("${CMAKE_Fortran_COMPILER_ID}" EQUAL "GNU")
    set(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} -fsanitize=address ")
  endif()
  if ("${CMAKE_C_COMPILER_ID}" EQUAL "GNU")
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} " -fsanitize=address ")
  endif()
  if ("${CMAKE_CXX_COMPILER_ID}" EQUAL "GNU")
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -fsanitize=address ")
  endif()
endif ()
# set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} " -fsanitize=address -DEP_PIO ")

if ( EP_PIO )
  set(C_HEADERS ${C_HEADERS}
    dependency/piotools/pioParallelUtils.h
    )
  set(C_SOURCES ${C_SOURCES}
    dependency/piotools/pioParallelUtils.c
    )
  set(CXX_SOURCES ${CXX_SOURCES}
    dependency/piotools/pioInterface.cpp
    )
  set(CXX_HEADERS ${CXX_HEADERS}
    dependency/piotools/pioInterface.hpp
    dependency/piotools/pio.hpp
    )
  set(F90_SOURCES ${F90_SOURCES}
    dependency/piotools/pio_interface.F90
    )
  
  set(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} -DEP_PIO=1 ")
  
endif ()

add_executable(mygrad ${F90_SOURCES} ${CXX_SOURCES} ${C_SOURCES} ${CXX_HEADERS} ${C_HEADERS})
set_property(TARGET mygrad PROPERTY LINKER_LANGUAGE Fortran)

if ( EP_GZ )
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -DENABLE_GZ=1 ")
  set_target_properties(mygrad PROPERTIES LINK_FLAGS " -lz ")
endif ()

if ( EP_MPI )
  find_package(MPI REQUIRED)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  target_link_libraries(mygrad ${MPI_Fortran_LIBRARIES})
  
  if(MPI_Fortran_COMPILE_FLAGS)
    set_target_properties(mygrad PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
  endif()

  if(MPI_Fortran_LINK_FLAGS)
    set_target_properties(mygrad PROPERTIES
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
  endif()


  # PIO mpi flag
  set(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} -DENABLE_MPI=1 ")
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -DENABLE_MPI=1 ")
  set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -DENABLE_MPI=1 ")
endif()
