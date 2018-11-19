cmake_minimum_required(VERSION 3.7)
project(sph3d-reworked)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -Wextra")

add_executable(sph3d main.cpp Cell.h Common.h Grid.h Params.h Particle.h ParticlesState.h Point.h Cell.cpp utils.h Grid.cpp ParticlesState.cpp Particle.cpp Common.cpp utils.cpp Numbered.cpp Numbered.h)