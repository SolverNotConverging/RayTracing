cmake_minimum_required(VERSION 4.2)
project(RayTracing)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
find_package(Matplot++ REQUIRED)

add_executable(RayTracing main.cpp)

target_link_libraries(RayTracing
        PRIVATE
        Eigen3::Eigen
        Matplot++::matplot
)



