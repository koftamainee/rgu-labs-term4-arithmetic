cmake_minimum_required(VERSION 3.25)

add_executable(viewer viewer.cpp)
target_include_directories(viewer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
