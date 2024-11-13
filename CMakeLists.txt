cmake_minimum_required(VERSION 3.15...3.27)

project(teqpflsh)
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

set(CMAKE_CXX_STANDARD 20)

set(BUILD_SHARED_LIBS TRUE "Use shared libraries for GEOS, due to the LGPL v2.1 license")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/externals/geos EXCLUDE_FROM_ALL)
set(BUILD_SHARED_LIBS FALSE "Everything else is statically linked")
    
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/externals/nanobind)

set(TEQP_NO_TESTS TRUE CACHE BOOL "No test targets for teqp (quite slow to compile)")
set(TEQP_NO_PYTHON TRUE CACHE BOOL "No pybind11 targets for teqp (quite slow to compile)")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp)

set(NANOFLANN_BUILD_EXAMPLES FALSE CACHE BOOL "No nanoflann examples")
set(NANOFLANN_BUILD_TESTS FALSE CACHE BOOL "No nanoflann tests")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/externals/nanoflann)

set(NAME teqpflsh)

# PYTHON STUFF
nanobind_add_module(_${NAME}_impl ${CMAKE_CURRENT_SOURCE_DIR}/interface/nanobind_interface.cpp)
target_include_directories(_${NAME}_impl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(_${NAME}_impl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp/externals/Eigen)
target_link_libraries(_${NAME}_impl PUBLIC geos PUBLIC teqpcpp PRIVATE nanoflann)

nanobind_add_stub(
  _${NAME}_impl_stub
  MODULE _${NAME}_impl
  OUTPUT "${CMAKE_SOURCE_DIR}/src/teqpflsh/_${NAME}_impl.pyi"
  PYTHON_PATH $<TARGET_FILE_DIR:_${NAME}_impl>
  DEPENDS _${NAME}_impl
  MARKER_FILE py.typed
  VERBOSE
)
# Install directive for scikit-build-core for the stubs
install(FILES "${CMAKE_SOURCE_DIR}/src/teqpflsh/_${NAME}_impl.pyi" DESTINATION ${NAME})
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/py.typed" DESTINATION ${NAME})

if (XCODE_DEBUG_PYTHON)
#See https://stackoverflow.com/a/72078672 and https://cmake.org/cmake/help/latest/manual/cmake-properties.7.html#properties-on-targets
set_target_properties(_${NAME}_impl PROPERTIES
                      XCODE_GENERATE_SCHEME TRUE
                      XCODE_SCHEME_WORKING_DIRECTORY "/Users/ihb/Documents/PapersInProgress/Working/(J)phasedet pure/figs"
                      XCODE_SCHEME_ARGUMENTS teqpflsh_Flasher.py
                      XCODE_SCHEME_EXECUTABLE /Users/ihb/mambaforge/envs/py310/bin/python)
endif()

# Install directive for scikit-build-core
install(TARGETS _${NAME}_impl LIBRARY DESTINATION ${NAME})
install(TARGETS geos LIBRARY DESTINATION ${NAME})

# Alter the output location for editable builds so that GEOS 
# is loaded properly
# TODO: do this in a better way; experiment-based
if (SKBUILD_STATE EQUAL "editable")
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set_target_properties( geos
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "$<TARGET_FILE_DIR:_${NAME}_impl>"
        LIBRARY_OUTPUT_DIRECTORY "$<TARGET_FILE_DIR:_${NAME}_impl>"
        RUNTIME_OUTPUT_DIRECTORY "$<TARGET_FILE_DIR:_${NAME}_impl>"
    )
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# Single-source the version, either from scikit, or from parsing the pyproject.toml
if (SKBUILD)
    add_definitions("-DTEQPFLSHVERSION=\"${SKBUILD_PROJECT_VERSION_FULL}\"")
else()
    file(READ "pyproject.toml" TOML_CONTENT)
    set(REG "version = \"([0-9]+\\.[0-9]+\\.[0-9]+)\"")
    string(REGEX MATCH "${REG}" VERSION_MATCH "${TOML_CONTENT}")
    if (NOT VERSION_MATCH)
        message(FATAL_ERROR "Can't parse the version")
    else()
        string(REGEX REPLACE "${REG}" "\\1" PROJECT_VERSION_FULL "${VERSION_MATCH}")
        message(STATUS "Version: ${PROJECT_VERSION_FULL}")
        add_definitions("-DTEQPFLSHVERSION=\"${PROJECT_VERSION_FULL}\"")
    endif()
endif()

if (TEQPFLSH_BUILD_TESTS)

    add_executable(geos_test "${CMAKE_CURRENT_SOURCE_DIR}/src/geos_test.cpp")
    target_include_directories(geos_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/externals/geos/include)
    target_link_libraries(geos_test PUBLIC geos)
    
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/externals/Catch2")
    add_executable(quadtree "${CMAKE_CURRENT_SOURCE_DIR}/src/quadtree.cpp")
    target_include_directories(quadtree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(quadtree PUBLIC geos PRIVATE Catch2WithMain PRIVATE teqpcpp PRIVATE nanoflann)

    add_executable(quadtree_orig "${CMAKE_CURRENT_SOURCE_DIR}/src/quadtree_orig.cpp")
    target_link_libraries(quadtree_orig PUBLIC geos PRIVATE Catch2WithMain)

    add_executable(geos_bench "${CMAKE_CURRENT_SOURCE_DIR}/src/geos_bench.cpp")
    target_link_libraries(geos_bench PUBLIC geos PRIVATE Catch2WithMain)
endif()
