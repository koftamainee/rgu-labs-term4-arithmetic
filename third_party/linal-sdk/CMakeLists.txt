add_library(linal STATIC
        src/core/bigfloat/arithmetics.cpp
        src/core/bigfloat/bigfloat.cpp
        src/core/bigfloat/exponential.cpp
        src/core/bigfloat/helpers.cpp
        src/core/bigfloat/io.cpp
        src/core/bigfloat/logarithmic.cpp
        src/core/bigfloat/logical.cpp
        src/core/bigfloat/trigonometry.cpp
        src/core/bigint/advanced.cpp
        src/core/bigint/arithmetics.cpp
        src/core/bigint/bitwise.cpp
        src/core/bigint/bigint.cpp
        src/core/bigint/division.cpp
        src/core/bigint/helpers.cpp
        src/core/bigint/io.cpp
        src/core/bigint/logical.cpp
        src/core/bigint/multiplication.cpp
        src/core/bigint/unsafe.cpp
        #        src/math/line_2d.cpp
        #        src/math/line_nd.cpp
        src/math/matrix.cpp
        src/math/vector.cpp
)


target_include_directories(linal PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include/core
        ${CMAKE_CURRENT_SOURCE_DIR}/include/math
)