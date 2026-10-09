# Use the project's existing environment for ordinary CMake builds too.
if (NOT DEFINED Python_EXECUTABLE AND EXISTS "${PROJECT_SOURCE_DIR}/.venv/Scripts/python.exe")
    set(Python_EXECUTABLE "${PROJECT_SOURCE_DIR}/.venv/Scripts/python.exe" CACHE FILEPATH "Python interpreter")
endif ()
find_package(Python 3.14 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
        OUTPUT_VARIABLE _pybind11_dir OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY)
file(TO_CMAKE_PATH "${_pybind11_dir}" _pybind11_dir)
find_package(pybind11 3.0 CONFIG REQUIRED HINTS "${_pybind11_dir}")
pybind11_add_module(_core MODULE NO_EXTRAS
        bindings/module.cpp bindings/geometry.cpp bindings/antennas.cpp
        bindings/solver.cpp bindings/results.cpp)
target_link_libraries(_core PRIVATE RayTracing::IO)
target_compile_definitions(_core PRIVATE RAYTRACING_VERSION="${PROJECT_VERSION}")

install(TARGETS _core RUNTIME_DEPENDENCY_SET raytracing_runtime
        LIBRARY DESTINATION raytracing RUNTIME DESTINATION raytracing)
install(FILES "${PROJECT_SOURCE_DIR}/antenna_patterns/simple_patch.ffs"
        DESTINATION raytracing/data)

# Ship the numerical/persistence runtime. Rendering dependencies are Python packages.
if (WIN32)
    install(RUNTIME_DEPENDENCY_SET raytracing_runtime
            DIRECTORIES "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin"
            PRE_EXCLUDE_REGEXES "api-ms-.*" "ext-ms-.*" "python[0-9]+(_d)?\\.dll"
            POST_EXCLUDE_REGEXES ".*[/\\\\][Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\\\].*"
            RUNTIME DESTINATION raytracing/.libs)
    # Include the license notices for the native dependencies supplied by vcpkg.
    foreach (_dependency IN ITEMS eigen3 nlohmann-json hdf5 zlib)
        install(FILES "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/${_dependency}/copyright"
                DESTINATION "raytracing/licenses/${_dependency}")
    endforeach ()
endif ()
