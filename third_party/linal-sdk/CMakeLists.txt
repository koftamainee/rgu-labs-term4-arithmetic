add_library(linal STATIC
        #        src/math/line_2d.cpp
        #        src/math/line_nd.cpp
        src/math/matrix.cpp
        src/math/vector.cpp
)

include(FetchContent)
FetchContent_Declare(
        bigmath
        GIT_REPOSITORY https://github.com/koftamainee/bigmath.git
        GIT_TAG master
)
set(BIGMATH_BUILD_TESTS   OFF CACHE BOOL "" FORCE)
set(BIGMATH_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(bigmath)

target_link_libraries(linal PUBLIC bigmath++)



target_include_directories(linal PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include/core
        ${CMAKE_CURRENT_SOURCE_DIR}/include/math
)