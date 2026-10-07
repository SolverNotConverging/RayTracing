cmake_minimum_required(VERSION 4.2)
project(RayTracing)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
find_package(VTK REQUIRED COMPONENTS
        CommonCore
        CommonDataModel
        InteractionStyle
        RenderingAnnotation
        RenderingCore
        RenderingFreeType
        RenderingOpenGL2
)

add_executable(RayTracing src/main.cpp
        include/EMRay.hpp
        src/EMRay.cpp
        src/Visualization.cpp
        include/Visualization.hpp)

target_include_directories(RayTracing
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(RayTracing
        PRIVATE
        Eigen3::Eigen
        ${VTK_LIBRARIES}
)

vtk_module_autoinit(TARGETS RayTracing MODULES ${VTK_LIBRARIES})
