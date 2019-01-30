cmake_minimum_required(VERSION 3.7)
project(sph3d-reworked)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wno-long-long -pedantic -Wextra")

add_executable(sph3d
                main.cpp

                Cell.h
                Common.h
                Grid.h
                Params.h
                Particle.h
        Point.h
                Numbered.h
                utils.h

                Cell.cpp
                Grid.cpp
        Particle.cpp
                Common.cpp
                utils.cpp
                Numbered.cpp
        CompareFiles.cpp CompareFiles.h Solver.h Solver.cpp InitStates.cpp InitStates.h)