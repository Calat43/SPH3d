cmake_minimum_required(VERSION 3.7)
project(sph3d-reworked)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wno-long-long -pedantic -Wextra")

add_executable(sph3d
                main.cpp

                Params.h
                Point.h
                Cell.h
                Common.h
                Grid.h
                Particle.h
                Numbered.h
                utils.h
                CompareFiles.h
                Solver.h
                InitStates.h

                Cell.cpp
                Common.cpp
                Grid.cpp
                Particle.cpp
                Numbered.cpp
                utils.cpp
                CompareFiles.cpp
                Solver.cpp
                InitStates.cpp SodTube3d.cpp SodTube3d.h Dusty_shock.cpp Dusty_shock.h Point.cpp NonLinear.cpp NonLinear.h DustyShock1d.cpp DustyShock1d.h)