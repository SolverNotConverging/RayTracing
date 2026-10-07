cmake_minimum_required(VERSION 4.2)
project(RayTracing)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(VTK REQUIRED COMPONENTS
        CommonCore
        CommonDataModel
        FiltersSources
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
        src/TraceConfig.cpp
        include/TraceConfig.hpp
        include/Visualization.hpp)

target_include_directories(RayTracing
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(RayTracing
        PRIVATE
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
        ${VTK_LIBRARIES}
)

vtk_module_autoinit(TARGETS RayTracing MODULES ${VTK_LIBRARIES})

# Default to the editable source file even when the IDE runs in the build folder.
target_compile_definitions(RayTracing PRIVATE
        RAYTRACING_DEFAULT_CONFIG="${CMAKE_CURRENT_SOURCE_DIR}/trace_options.json")
