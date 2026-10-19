cmake_minimum_required(VERSION 3.20)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(dna_sequence_similarity LANGUAGES CXX)

file(GLOB_RECURSE CPP_SOURCES CONFIGURE_DEPENDS
	"${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(debug src/main.cpp)
target_sources(debug PRIVATE ${CPP_SOURCES})
target_compile_features(debug PRIVATE cxx_std_23)
target_compile_options(debug PRIVATE
	-Wall
	-Wextra
	-march=native
	-Wshadow
	-ggdb
)

add_executable(run src/main.cpp)
target_sources(run PRIVATE ${CPP_SOURCES})
target_compile_features(run PRIVATE cxx_std_23)
target_compile_options(run PRIVATE
	-O2
	-march=native
)

# cmake --build build --target=format
file(GLOB_RECURSE ALL_CXX_SOURCE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

add_custom_target(format
    COMMAND clang-format -i ${ALL_CXX_SOURCE_FILES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Formatting all source files with clang-format"
)
