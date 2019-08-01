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
  "../sources/base"
)

set(Sph3dFilesBase
  "../sources/base/Point.h"
  "../sources/base/Point.cpp"
  "../sources/base/Cell.h"
  "../sources/base/Cell.cpp"
  "../sources/base/MathUtils.h"
  "../sources/base/MathUtils.cpp"
  "../sources/base/Grid.h"
  "../sources/base/Grid.cpp"
  "../sources/base/Particle.h"
  "../sources/base/Particle.cpp"
  "../sources/base/Numbered.h"
  "../sources/base/Numbered.cpp"
  "../sources/base/Viscosity.h"
  "../sources/base/Viscosity.cpp"
)

set(Sph3dFilesSolvers
  "../sources/CompareFiles.h"
  "../sources/CompareFiles.cpp"
  "../sources/SodTube1d.h"
  "../sources/SodTube1d.cpp"
  "../sources/InitStates.h"
  "../sources/InitStates.cpp"
  "../sources/SodTube3d.h"
  "../sources/SodTube3d.cpp"
  "../sources/DustyShock.h"
  "../sources/DustyShock.cpp"
  "../sources/NonLinear.h"
  "../sources/NonLinear.cpp"
  "../sources/DustyShock1d.h"
  "../sources/DustyShock1d.cpp"
  "../sources/BallInVacuum.h"
  "../sources/BallInVacuum.cpp"
  "../sources/SPHSolver.h"
  "../sources/SPHSolver.cpp"
)

set(Sph3dFilesSettings
  "../sources/Params.h"
)


source_group("base" FILES ${Sph3dFilesBase})
source_group("solvers" FILES ${Sph3dFilesSolvers})
source_group("settings" FILES ${Sph3dFilesSettings})

set(Sph3dLib
  ${Sph3dFilesBase}
  ${Sph3dFilesSolvers}
  ${Sph3dFilesSettings}
)

add_library(Sph3d STATIC ${Sph3dLib})

add_executable(Sph3dTestsRunner "../sources/SPHTestsRunner.cpp")
target_link_libraries(Sph3dTestsRunner Sph3d)
if (NOT DEFINED MSVC)
  target_link_libraries(Sph3dTestsRunner m)  # link libm math library
endif()

add_executable(Sph3dSolverRunner "../sources/SPHSolverRunner.cpp")
target_link_libraries(Sph3dSolverRunner Sph3d)
if (NOT DEFINED MSVC)
  target_link_libraries(Sph3dSolverRunner m)  # link libm math library
endif()


