cmake_minimum_required(VERSION 3.3)

project(satid C)

set(CMAKE_C_FLAGS "-std=c99 -Wall -Wextra")

find_package(PkgConfig)
pkg_check_modules(LIBPREDICT predict)
include_directories(${LIBPREDICT_INCLUDE_DIRS})

add_executable(satid main.c)
target_link_libraries(satid m ${LIBPREDICT_LIBRARIES})
