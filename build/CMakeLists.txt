cmake_minimum_required(VERSION 2.8.1)
project(sph3d-reworked)

set(CMAKE_CXX_STANDARD 11)

# define platform
if (MSVC)
  set(PLATFORM "VS2017_Win64")
else()
  set(PLATFORM "Linux_x64")
endif()

# define compiler options
if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  string(REPLACE "INCREMENTAL:YES" "INCREMENTAL:NO" replacementFlags ${CMAKE_EXE_LINKER_FLAGS_DEBUG}) 
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/INCREMENTAL:NO ${replacementFlags}" )
  string(REPLACE "INCREMENTAL:YES" "INCREMENTAL:NO" replacementFlags ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}) 
  set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "/INCREMENTAL:NO ${replacementFlags}" )

  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  add_definitions(-D_WINDOWS_SOURCE)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-long-long -Wextra -Wall -pedantic")
endif()

# define directories (put all configurations into single ./bin dir)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY                "../bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG          "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE        "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG          "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE        "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG          "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE        "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")


file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

include_directories(
  "../sources"
)

set(Sph3dFiles
  "../sources/Params.h"
  "../sources/Point.h"
  "../sources/Point.cpp"
  "../sources/Cell.h"
  "../sources/Cell.cpp"
  "../sources/Common.h"
  "../sources/Common.cpp"
  "../sources/Grid.h"
  "../sources/Grid.cpp"
  "../sources/Particle.h"
  "../sources/Particle.cpp"
  "../sources/Numbered.h"
  "../sources/Numbered.cpp"
  "../sources/utils.h"
  "../sources/utils.cpp"
  "../sources/CompareFiles.h"
  "../sources/CompareFiles.cpp"
  "../sources/Solver.h"
  "../sources/Solver.cpp"
  "../sources/InitStates.h"
  "../sources/InitStates.cpp"
  "../sources/SodTube3d.h"
  "../sources/SodTube3d.cpp"
  "../sources/Dusty_shock.h"
  "../sources/Dusty_shock.cpp"
  "../sources/NonLinear.h"
  "../sources/NonLinear.cpp"
  "../sources/DustyShock1d.h"
  "../sources/DustyShock1d.cpp"
)

source_group("sources" FILES ${Sph3dFiles})

set(Sph3dLib
  ${Sph3dFiles}
)

add_library(Sph3d STATIC ${Sph3dLib})

add_executable(Sph3dRunner "../sources/main.cpp")
target_link_libraries(Sph3dRunner Sph3d)
if (NOT DEFINED MSVC)
  target_link_libraries(Sph3dRunner m)  # link libm math library
endif()
