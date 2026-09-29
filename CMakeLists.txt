cmake_minimum_required(VERSION 3.20)
project(sph_sim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -Wall -Wextra")

add_executable(sph_sim
    src/main.cpp
    src/Particles.cpp
)

if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    set(OpenMP_CXX_FLAGS
        "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY
        "${HOMEBREW_LIBOMP_PREFIX}/lib/libomp.dylib")
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(sph_sim PRIVATE OpenMP::OpenMP_CXX)

find_package(xsimd REQUIRED CONFIG)
target_link_libraries(sph_sim PRIVATE xsimd)