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
        src/Surface.cpp
        src/SequenceEvaluator.cpp
        src/PathRefinement.cpp
        include/PathRefinement.hpp
        include/SequenceEvaluator.hpp
        src/Visualization.cpp
        src/TraceConfig.cpp
        include/TraceConfig.hpp
        include/Visualization.hpp
        include/Surface.hpp)

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

include(CTest)
if (BUILD_TESTING)
    add_executable(SurfaceTests tests/SurfaceTests.cpp src/Surface.cpp src/EMRay.cpp)
    target_include_directories(SurfaceTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(SurfaceTests PRIVATE Eigen3::Eigen)
    add_test(NAME SurfaceTests COMMAND SurfaceTests)

    add_executable(SequenceEvaluatorTests tests/SequenceEvaluatorTests.cpp
            src/SequenceEvaluator.cpp src/Surface.cpp)
    target_include_directories(SequenceEvaluatorTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(SequenceEvaluatorTests PRIVATE Eigen3::Eigen)
    add_test(NAME SequenceEvaluatorTests COMMAND SequenceEvaluatorTests)
    add_executable(PathRefinementTests tests/PathRefinementTests.cpp
            src/PathRefinement.cpp src/SequenceEvaluator.cpp src/Surface.cpp src/EMRay.cpp)
    target_include_directories(PathRefinementTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(PathRefinementTests PRIVATE Eigen3::Eigen)
    add_test(NAME PathRefinementTests COMMAND PathRefinementTests)
endif ()
