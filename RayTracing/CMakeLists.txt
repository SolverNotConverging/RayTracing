find_package(Eigen3 CONFIG REQUIRED)
add_library(RayTracingCore STATIC
        src/Surface.cpp src/EMRay.cpp src/SequenceEvaluator.cpp src/PathRefinement.cpp
        src/Spreading.cpp src/FieldReconstruction.cpp src/ImpulseResponse.cpp
        src/Antenna.cpp src/Solver.cpp)
add_library(RayTracing::RayTracing ALIAS RayTracingCore)
target_include_directories(RayTracingCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include/RayTracing)
target_link_libraries(RayTracingCore PUBLIC Eigen3::Eigen)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(hdf5 CONFIG REQUIRED)
add_library(RayTracingIO STATIC src/ResultIO.cpp)
target_compile_definitions(RayTracingIO PRIVATE RAYTRACING_VERSION="${PROJECT_VERSION}")
add_library(RayTracing::IO ALIAS RayTracingIO)
target_link_libraries(RayTracingIO PUBLIC RayTracingCore PRIVATE nlohmann_json::nlohmann_json hdf5::hdf5-shared)
