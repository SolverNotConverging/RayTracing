find_package(Python 3.14 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
        OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
find_package(pybind11 3.0 CONFIG REQUIRED)

pybind11_add_module(_core MODULE
        bindings/module.cpp bindings/common.cpp bindings/point_ray.cpp bindings/beam.cpp)
target_link_libraries(_core PRIVATE RayTracing::RayTracing)
target_compile_definitions(_core PRIVATE RAYTRACING_VERSION="${PROJECT_VERSION}")
install(TARGETS _core LIBRARY DESTINATION raytracing RUNTIME DESTINATION raytracing)

# Installation can strip the module and invalidate its ad-hoc signature on macOS arm64.
if (APPLE)
    install(CODE [[
        execute_process(COMMAND /usr/bin/codesign --force --sign -
                "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/raytracing/$<TARGET_FILE_NAME:_core>"
                COMMAND_ERROR_IS_FATAL ANY)
    ]])
endif ()
