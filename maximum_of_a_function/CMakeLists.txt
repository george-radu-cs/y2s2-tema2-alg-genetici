cmake_minimum_required(VERSION 3.21)
project (Tema_2___alg_genetici_maximum_of_a_function)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

add_library(
        genetic_algorithm
        src/RandomGenerator.h
        src/Chromosome.h
        src/Chromosome.cpp
        src/GeneticAlgorithm.h
        src/GeneticAlgorithm.cpp
)

add_executable(Tema_2___alg_genetici_maximum_of_a_function src/main.cpp)

target_link_libraries(Tema_2___alg_genetici_maximum_of_a_function PRIVATE genetic_algorithm)
