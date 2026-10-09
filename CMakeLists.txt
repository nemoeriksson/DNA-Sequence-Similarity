cmake_minimum_required(VERSION 3.20)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(run LANGUAGES CXX)

file(GLOB_RECURSE CPP_SOURCES CONFIGURE_DEPENDS
	"${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(run src/main.cpp)

target_sources(run
	PRIVATE
	${CPP_SOURCES}
)

target_compile_features(run PRIVATE cxx_std_23)

target_compile_options(run PRIVATE
	-Wall
	-Wextra
	-march=native
	-Wshadow
	-ggdb
)

