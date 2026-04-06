cmake_minimum_required(VERSION 3.25)
project(rgu-labs-term4-arithmetic LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")

find_program(CCACHE ccache)
if(CCACHE)
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
endif()

include_directories(include)

include(FetchContent)
FetchContent_Declare(
        chalkboard
        GIT_REPOSITORY https://github.com/koftamainee/chalkboard.git
        GIT_TAG        master
)
FetchContent_MakeAvailable(chalkboard)
set(CHALKBOARD_BUILD_VIEWER ON CACHE BOOL "" FORCE)

file(GLOB TASK_SOURCES "src/*.cpp")
foreach(src_file ${TASK_SOURCES})
    get_filename_component(task_name ${src_file} NAME_WE)
    if (task_name STREQUAL "631")
        continue()
    endif ()
    add_executable(${task_name} ${src_file})
    target_link_libraries(${task_name} gmp gmpxx chalkboard::chalkboard)
    target_precompile_headers(${task_name} PRIVATE include/pch.hpp)
endforeach()

include(FetchContent)

FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG 3.4
)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

FetchContent_Declare(
        imgui
        GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG v1.91.6
)
FetchContent_MakeAvailable(imgui)

add_library(imgui_lib STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui_lib PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui_lib PUBLIC glfw)

find_package(OpenGL REQUIRED)

add_executable(631 src/631.cpp)
target_include_directories(631 PRIVATE include)
target_link_libraries(631 PRIVATE imgui_lib OpenGL::GL)