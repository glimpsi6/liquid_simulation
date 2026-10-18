cmake_minimum_required(VERSION 3.20)
project(sph_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ENABLE_ASAN   "Build with AddressSanitizer + UBSan" OFF)
option(ENABLE_NATIVE "Use native CPU tuning flags (non-portable binary)" OFF)

add_executable(sph_sim
    src/main.cpp
    src/Particles.cpp
)
target_include_directories(sph_sim PRIVATE include)

target_compile_options(sph_sim PRIVATE -Wall -Wextra
                                       $<$<CONFIG:Release>:-g>)

if(ENABLE_NATIVE AND NOT ENABLE_ASAN)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    target_compile_options(sph_sim PRIVATE $<$<CONFIG:Release>:-mcpu=native>)
  else()
    target_compile_options(sph_sim PRIVATE $<$<CONFIG:Release>:-march=native>)
  endif()
endif()

# Санитайзеры
if(ENABLE_ASAN)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(sph_sim PRIVATE
        -fsanitize=address,undefined -fno-omit-frame-pointer -g)
    target_link_options(sph_sim PRIVATE -fsanitize=address,undefined)
  else()
    message(WARNING "ENABLE_ASAN: компилятор не Clang/GCC, опция пропущена")
  endif()
endif()

# OpenMP
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  find_program(BREW_EXECUTABLE brew)
  if(BREW_EXECUTABLE)
    execute_process(
        COMMAND ${BREW_EXECUTABLE} --prefix libomp
        OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
  endif()
  if(HOMEBREW_LIBOMP_PREFIX)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${HOMEBREW_LIBOMP_PREFIX}/lib/libomp.dylib")
  else()
    message(WARNING "libomp не найден: выполните brew install libomp")
  endif()
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(sph_sim PRIVATE OpenMP::OpenMP_CXX)