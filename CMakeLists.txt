cmake_minimum_required(VERSION 3.10)

# Set the project name
project(Legacy)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable clang-tidy with cppcoreguideline warnings
find_program(CLANG_TIDY "clang-tidy")
if(CLANG_TIDY)
  set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY};-checks=cppcoreguidelines-*")
endif()

# Add the executable
add_subdirectory(naca)
