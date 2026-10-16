# Builds inside this project, or standalone against an installed package:
#   cmake -S examples/cpp -B build/examples -DCMAKE_PREFIX_PATH=<install prefix>
cmake_minimum_required(VERSION 4.2)
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(RayTracingExamples LANGUAGES CXX)
    find_package(RayTracing 1 CONFIG REQUIRED)
endif ()
foreach (example IN ITEMS multipath gaussian_beam)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE RayTracing::RayTracing)
    if (BUILD_TESTING)
        add_test(NAME example_${example} COMMAND ${example})
    endif ()
endforeach ()
