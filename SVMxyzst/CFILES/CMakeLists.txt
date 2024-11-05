cmake_minimum_required(VERSION 3.14)
project(SVMmagnetic)

set(CMAKE_CXX_STANDARD 14)

include_directories(/home/kirscher/CLionProjects/SVMmagnetic/eigen/Eigen)

add_library(eigen3 SHARED IMPORTED)
add_executable(SVMmagnetic Input.h Input.cpp BasisState.h BasisState.cpp MatrixElement.h MatrixElement.cpp Rand.h Rand.cpp Permutation.h Permutation.cpp
        SVM.h SVM.cpp Source.cpp Operators.h Operators.cpp CoordinatsTransformation.h CoordinatsTransformation.cpp)