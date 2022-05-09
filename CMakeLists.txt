######################################
# CMakeLists.txt for MCNP QT interface
######################################
cmake_minimum_required(VERSION 3.7.0)

project(mygrad LANGUAGES Fortran)

# Specify default build
set(CMAKE_BUILD_TYPE Release CACHE STRING "Buld type")

# Specify sources
set(F90_SOURCES
  src/clone_lib.f90
  src/define_kind.f90 
  src/fakemesh.f90 
  src/mesh_types.f90
  src/util.f90
  src/derivatives_common_template.f90
  )

set_source_files_properties(
  ${F90_SOURCES}
  PROPERTIES Fortran_PREPROCESS ON
)

add_executable(mygrad ${F90_SOURCES} )
