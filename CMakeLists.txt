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

# find_package(Threads REQUIRED)
# target_link_libraries(sph_sim PRIVATE Threads::Threads)