find_package(Eigen3 CONFIG REQUIRED)
add_library(RayTracingCore STATIC
        src/Surface.cpp src/Material.cpp src/Tracing.cpp src/SequenceEvaluator.cpp src/PathRefinement.cpp
        src/Spreading.cpp src/FieldReconstruction.cpp src/ImpulseResponse.cpp
        src/Antenna.cpp src/Solver.cpp src/BeamSolver.cpp)
add_library(RayTracing::Core ALIAS RayTracingCore)
set_target_properties(RayTracingCore PROPERTIES EXPORT_NAME Core POSITION_INDEPENDENT_CODE ON)
target_compile_features(RayTracingCore PUBLIC cxx_std_20)
target_include_directories(RayTracingCore PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(RayTracingCore PUBLIC Eigen3::Eigen)
install(TARGETS RayTracingCore EXPORT RayTracingCoreTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
install(DIRECTORY include/RayTracing DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT Development PATTERN ResultIO.hpp EXCLUDE)
install(EXPORT RayTracingCoreTargets NAMESPACE RayTracing::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RayTracing COMPONENT Development)

if (RAYTRACING_BUILD_IO)
    enable_language(C)
    find_package(nlohmann_json CONFIG REQUIRED)
    set(HDF5_USE_STATIC_LIBRARIES OFF)
    find_package(HDF5 MODULE REQUIRED COMPONENTS C)
    add_library(RayTracingIO STATIC src/ResultIO.cpp)
    add_library(RayTracing::IO ALIAS RayTracingIO)
    set_target_properties(RayTracingIO PROPERTIES EXPORT_NAME IO POSITION_INDEPENDENT_CODE ON)
    target_compile_definitions(RayTracingIO PRIVATE RAYTRACING_VERSION="${PROJECT_VERSION}")
    target_link_libraries(RayTracingIO PUBLIC RayTracing::Core
            PRIVATE nlohmann_json::nlohmann_json HDF5::HDF5)
    install(TARGETS RayTracingIO EXPORT RayTracingIOTargets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
    install(FILES include/RayTracing/ResultIO.hpp
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/RayTracing COMPONENT Development)
    install(EXPORT RayTracingIOTargets NAMESPACE RayTracing::
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RayTracing COMPONENT Development)
endif ()
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/RayTracingConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/RayTracingConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RayTracing)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/RayTracingConfigVersion.cmake
        VERSION ${PROJECT_VERSION} COMPATIBILITY SameMinorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/RayTracingConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/RayTracingConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RayTracing COMPONENT Development)
