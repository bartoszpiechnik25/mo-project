cmake_minimum_required(VERSION 3.23)
project(ProjektLab11)

set(CMAKE_CXX_STANDARD 20)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")


add_executable(ProjektLab11 main.cpp)

#target_sources(ProjektLab11 PRIVATE "${SRC_DIR}/stale.cpp")
target_include_directories(ProjektLab11 PRIVATE ${INC_DIR})
