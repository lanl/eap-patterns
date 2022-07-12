######################################
# CMakeLists.txt for EAP Patterns
######################################
cmake_minimum_required(VERSION 3.7.0)
project(mygrad LANGUAGES Fortran CXX)

option(EP_PIO "Build targets with PIO integration;" "ON")

# unimplemented options
option(EP_HDF5 "Build targets with HDF5 integration; not implemented yet" "OFF")
option(EP_MPI "Build targets with MPI integration; not implemented yet" "OFF")
option(EP_OPENMP "Build targets with OpenMP integration; not implemented yet" "OFF")

# Specify default build
set(CMAKE_BUILD_TYPE Release CACHE STRING "Buld type")

# Specify sources
set(F90_SOURCES
  src/clone_lib.f90
  src/define_kind.f90 
  src/fakemesh.f90 
  src/mem_release.f90
  src/mesh_types.f90
  src/util.f90
  src/main.F90
  src/derivatives_common_template.f90
  )

set_source_files_properties(
  ${F90_SOURCES}
  PROPERTIES Fortran_PREPROCESS ON
)

if ( EP_PIO )
  set(CXX_SOURCES ${CXX_SOURCES}
    dependency/piotools/pioInterface.cpp
    )
  set(CXX_HEADERS ${CXX_HEADERS}
    dependency/piotools/pioInterface.hpp
    dependency/piotools/pio.hpp
    )
  set(CXX_SOURCES ${CXX_SOURCES}
    dependency/piotools/pio_interface.F90
    )
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} -DEP_PIO)
endif ()

add_executable(mygrad ${F90_SOURCES} ${CXX_SOURCES} ${CXX_HEADERS})
