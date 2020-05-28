cmake_minimum_required(VERSION 3.17)

set(CMAKE_C_FLAG "${CMAKE_C_FLAGS} -g -O0 -fmax-errors=5 -Wfatal-errors")
set(CMAKE_CXX_FLAG "${CMAKE_C_FLAGS} -g -O0 -fmax-errors=5 -Wfatal-errors -enable-threads")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(HOGE LANGUAGES C CXX)

find_package(FFTW3 REQUIRED)

find_library (
    FFTW3_LIB
    NAMES fftw3 libfftw3 # what to look for
    HINTS ${FFTW3_LIBRARY_DIRS} # where to look
    NO_DEFAULT_PATH # do not search system default paths
)

find_library (
    FFTW3_THREAD_LIB
    NAMES fftw3_threads libfftw3_threads # what to look for
    HINTS ${FFTW3_LIBRARY_DIRS} # where to look
    NO_DEFAULT_PATH # do not search system default paths
)

set(EXEC_NAMES isolated test)

foreach(EXEC_NAME ${EXEC_NAMES})

    add_executable(${EXEC_NAME} "${EXEC_NAME}")
    target_include_directories(${EXEC_NAME} PRIVATE
        ${FFTW3_INCLUDE_DIRS}
        external/grid/include
        external/unit/include)
    target_link_libraries(${EXEC_NAME} ${FFTW3_LIB} ${FFTW3_THREAD_LIB} m)

endforeach()
