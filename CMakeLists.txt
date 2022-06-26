cmake_minimum_required(VERSION 3.10)

# set the project name
project(ParkPlacement)

# bring headers into project
file (STRINGS "routingkit_include.txt" _RK_INCLUDE)
file (STRINGS "routingkit_lib.txt" _RK_LIB)
include_directories(include)
include_directories(${_RK_INCLUDE})
# set sources
file(GLOB SOURCES "src/*.cpp")

# add the executable
add_executable(ParkPlacement ${SOURCES})
target_link_directories(ParkPlacement PUBLIC ${_RK_LIB})

# link library needed for the vertex parsing
target_link_libraries(ParkPlacement pthread)
target_link_libraries(ParkPlacement routingkit)
target_link_libraries(ParkPlacement boost_filesystem)
target_link_libraries(ParkPlacement boost_program_options)

# specify c++ standard
set_property(TARGET ParkPlacement PROPERTY CXX_STANDARD 17)