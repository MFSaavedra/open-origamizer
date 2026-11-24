# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

cmake_minimum_required(VERSION 3.28)
project(origamizer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ORIGAMIZER_BUILD_TESTS "Build the Catch2 unit tests" ON)
option(ORIGAMIZER_BUILD_PYTHON "Build the nanobind module origamizer._core" OFF)
option(ORIGAMIZER_SANITIZE "Instrument with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

# HiGHS (MIT): the LP and QP solver of the layout step (CLAUDE.md, Decisions), built
# from source at a pinned release as a static, position-independent library. It is
# added before the sanitizer flags, so it is not instrumented, and with
# warnings-as-errors off, so its warnings cannot fail our strict builds.
include(FetchContent)
set(_origamizer_werror ${CMAKE_COMPILE_WARNING_AS_ERROR})
set(CMAKE_COMPILE_WARNING_AS_ERROR OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
foreach(option BUILD_CXX_EXE BUILD_TESTING BUILD_EXAMPLES BUILD_SHARED_EXTRAS_LIB ZLIB)
  set(${option} OFF CACHE BOOL "" FORCE)
endforeach()
FetchContent_Declare(highs
  GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
  GIT_TAG v1.15.1
  GIT_SHALLOW TRUE
  EXCLUDE_FROM_ALL
  SYSTEM)
FetchContent_MakeAvailable(highs)
set(CMAKE_COMPILE_WARNING_AS_ERROR ${_origamizer_werror})

if(ORIGAMIZER_SANITIZE)
  if(ORIGAMIZER_BUILD_PYTHON)
    message(FATAL_ERROR "ORIGAMIZER_SANITIZE is for the C++ tests; build the Python module without it.")
  endif()
  add_compile_options(-fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

# Debug builds are deliberate here (sanitizer runs), so skip CGAL's warning about them.
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
find_package(CGAL 5.6 REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

set(ORIGAMIZER_WARNINGS -Wall -Wextra -Wpedantic)

# Core library: every algorithm lives here (CLAUDE.md, C++ stack).
add_library(origamizer_core STATIC
  src/check.cpp
  src/cut.cpp
  src/detail/quadratic_program.cpp
  src/layout.cpp
  src/mesh.cpp
  src/molecules.cpp
)
target_include_directories(origamizer_core PUBLIC ${PROJECT_SOURCE_DIR}/include
                                           PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(origamizer_core PUBLIC CGAL::CGAL Eigen3::Eigen PRIVATE highs::highs)
target_compile_options(origamizer_core PRIVATE ${ORIGAMIZER_WARNINGS})

if(ORIGAMIZER_BUILD_TESTS)
  find_package(Catch2 3 REQUIRED)
  enable_testing()
  add_executable(origamizer_tests
    tests/cpp/test_check.cpp
    tests/cpp/test_check_ndebug.cpp
    tests/cpp/test_cut.cpp
    tests/cpp/test_face_shape.cpp
    tests/cpp/test_layout.cpp
    tests/cpp/test_molecules.cpp
    tests/cpp/test_obj.cpp
    tests/cpp/test_quadratic_program.cpp
    tests/cpp/test_soup.cpp
    tests/cpp/test_toolchain.cpp
    tests/cpp/test_topology.cpp
  )
  target_link_libraries(origamizer_tests PRIVATE origamizer_core Catch2::Catch2WithMain)
  # Tests may include private headers under src/ (for example src/detail).
  target_include_directories(origamizer_tests PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_compile_definitions(origamizer_tests PRIVATE ORIGAMIZER_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  target_compile_options(origamizer_tests PRIVATE ${ORIGAMIZER_WARNINGS})
  include(Catch)
  catch_discover_tests(origamizer_tests)
endif()

if(ORIGAMIZER_BUILD_PYTHON)
  find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC src/bindings/module.cpp)
  target_link_libraries(_core PRIVATE origamizer_core)
  install(TARGETS _core LIBRARY DESTINATION origamizer)
endif()
