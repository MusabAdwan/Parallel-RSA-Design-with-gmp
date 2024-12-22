cmake_minimum_required(VERSION 3.10)

project(MyProject)

set(CMAKE_CXX_STANDARD 11)

# Specify the vcpkg toolchain file
set(CMAKE_TOOLCHAIN_FILE "C:/Users/Dell/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")

find_package(GMP REQUIRED)

add_executable(my_executable main.cpp)

target_link_libraries(my_executable PRIVATE gmp)