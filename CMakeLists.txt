cmake_minimum_required(VERSION 4.2)
project(RayTracing VERSION 1.0.1 LANGUAGES CXX)

option(RAYTRACING_BUILD_PYTHON "Build the Python extension module" OFF)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
find_package(Eigen3 CONFIG REQUIRED)

add_library(raytracing STATIC
        src/antenna.cpp src/beam.cpp src/field.cpp src/geometry.cpp src/material.cpp
        src/refinement.cpp src/sequence.cpp src/solver.cpp src/spreading.cpp src/tracing.cpp)
add_library(RayTracing::RayTracing ALIAS raytracing)
set_target_properties(raytracing PROPERTIES EXPORT_NAME RayTracing POSITION_INDEPENDENT_CODE ON)
target_compile_features(raytracing PUBLIC cxx_std_20)
target_include_directories(raytracing PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(raytracing PUBLIC Eigen3::Eigen)

if (RAYTRACING_BUILD_PYTHON)
    add_subdirectory(python)
else ()
    install(TARGETS raytracing EXPORT RayTracingTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(DIRECTORY include/raytracing DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(EXPORT RayTracingTargets NAMESPACE RayTracing::
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RayTracing)
    configure_package_config_file(cmake/RayTracingConfig.cmake.in RayTracingConfig.cmake
            INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RayTracing)
    write_basic_package_version_file(RayTracingConfigVersion.cmake COMPATIBILITY SameMajorVersion)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/RayTracingConfig.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/RayTracingConfigVersion.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RayTracing)
endif ()

if (PROJECT_IS_TOP_LEVEL AND NOT RAYTRACING_BUILD_PYTHON)
    include(CTest)
    if (BUILD_TESTING)
        add_subdirectory(tests/cpp)
        add_subdirectory(examples/cpp)
    endif ()
endif ()
