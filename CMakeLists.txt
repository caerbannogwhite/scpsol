cmake_minimum_required(VERSION 3.18)
project(scpsol LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(highs
    GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
    GIT_TAG v1.9.0
    GIT_SHALLOW TRUE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(highs)

add_executable(scpsol
    src/main.cpp
    src/solver.cpp
    src/reader.cpp
    src/preprocessor.cpp
    src/lp.cpp
    src/heuristics.cpp
    src/cuts.cpp
    src/balas.cpp
    src/bnb.cpp
    src/reliability.cpp
    src/lagrangian.cpp
    src/decomposition.cpp)
target_include_directories(scpsol PRIVATE src)
target_link_libraries(scpsol PRIVATE highs)

if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".js")
    target_link_options(scpsol PRIVATE
        -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=67108864
        -sSTACK_SIZE=1048576 -sNODERAWFS=1 -sEXIT_RUNTIME=1)
endif()
