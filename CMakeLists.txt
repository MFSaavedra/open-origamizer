# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

cmake_minimum_required(VERSION 3.26)
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
)
target_include_directories(origamizer_core PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(origamizer_core PUBLIC CGAL::CGAL Eigen3::Eigen)
target_compile_options(origamizer_core PRIVATE ${ORIGAMIZER_WARNINGS})
set_target_properties(origamizer_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(ORIGAMIZER_BUILD_TESTS)
  find_package(Catch2 3 REQUIRED)
  enable_testing()
  add_executable(origamizer_tests
    tests/cpp/test_check.cpp
    tests/cpp/test_check_ndebug.cpp
    tests/cpp/test_toolchain.cpp
  )
  target_link_libraries(origamizer_tests PRIVATE origamizer_core Catch2::Catch2WithMain)
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
