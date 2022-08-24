######################################
# CMakeLists.txt for EAP Patterns
######################################
cmake_minimum_required(VERSION 3.16.0)
project(mygrad LANGUAGES Fortran C)

option(EP_MPI "Build targets with MPI integration;" "ON")
option(EP_SANITIZE "Add address sanitizer flag;" "ON")
option(EP_GZ "Enable GZ file reading in PIO;" "ON")

# Specify default build
set(CMAKE_BUILD_TYPE Release CACHE STRING "Buld type")

# Specify sources
set(F90_SOURCES
  src/binReader.f90
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

set(C_SOURCES
  src/fileUtils.c
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
endif ()


add_executable(mygrad ${F90_SOURCES} ${C_SOURCES} )
set_property(TARGET mygrad PROPERTY LINKER_LANGUAGE Fortran)

# find_package(VTUNE REQUIRED)

include_directories(/usr/projects/hpcsoft/toss3/common/x86_64/oneapi/2021.2.0.2883/vtune/2021.2.0/include/intel64)
# link_directories(/usr/projects/hpcsoft/toss3/common/x86_64/oneapi/2021.2.0.2883/vtune/2021.2.0/lib64)
target_link_libraries(mygrad "-L/usr/projects/hpcsoft/toss3/common/x86_64/oneapi/2021.2.0.2883/vtune/2021.2.0/lib64 -littnotify")

if ( EP_GZ )
  find_package(ZLIB REQUIRED)
  set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -DENABLE_GZ=1 ")
  include_directories(${ZLIB_INCLUDE_DIRS})
  target_link_libraries(mygrad ${ZLIB_LIBRARIES})
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
  set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -DENABLE_MPI=1 ")
endif()


include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

if( supported )
    message(STATUS "IPO / LTO enabled")
    set_property(TARGET mygrad PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()

