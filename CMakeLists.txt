cmake_minimum_required(VERSION 4.2)
project(RayTracing)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
find_package(Matplot++ REQUIRED)

add_executable(RayTracing src/main.cpp
        include/crt/ray.hpp
        include/crt/geometry.hpp
        include/crt/propagation.hpp
        include/crt/visualization.hpp)

target_include_directories(RayTracing
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(RayTracing
        PRIVATE
        Eigen3::Eigen
        Matplot++::matplot
)



