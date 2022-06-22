cmake_minimum_required(VERSION 3.10)

# set the project name
project(ParkPlacement)

# bring headers into project
include_directories(include)

# set sources
file(GLOB SOURCES "src/*.cpp")

# add the executable
add_executable(ParkPlacement ${SOURCES})
target_link_directories(ParkPlacement PUBLIC "lib")

# link library needed for the vertex parsing
target_link_libraries(ParkPlacement pthread)
target_link_libraries(ParkPlacement routingkit)
target_link_libraries(ParkPlacement boost_filesystem)
target_link_libraries(ParkPlacement boost_program_options)

# specify c++ standard
set_property(TARGET ParkPlacement PROPERTY CXX_STANDARD 17)