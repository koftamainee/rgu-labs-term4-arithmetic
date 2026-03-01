cmake_minimum_required(VERSION 3.25)
project(rgu-labs-term4-arithmetic LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

include_directories(include)

add_subdirectory(third_party/linal-sdk)

file(GLOB TASK_SOURCES "src/*.cpp")
foreach (src_file ${TASK_SOURCES})
    get_filename_component(task_name ${src_file} NAME_WE)
    add_executable(${task_name} ${src_file})
    target_link_libraries(${task_name} linal)
endforeach ()
