cmake_minimum_required(VERSION 3.22)
project(ePIC++ LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O2")

# Option to enable/disable plotting
option(ENABLE_PLOTTING "Enable Matplotlib plotting" OFF)

# Source files
file(GLOB SOURCES "src/*.cpp" "cross_section/*.cpp")

# Find HDF5 with CXX component (required for all builds)
find_package(HDF5 COMPONENTS CXX REQUIRED)
if(NOT HDF5_FOUND)
    message(FATAL_ERROR "HDF5 not found. Please install HDF5.")
endif()

# Find Threads (required for all builds)
find_package(Threads REQUIRED)

# Find FFTW3 (required for all builds)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)
if(NOT FFTW_FOUND)
    message(FATAL_ERROR "FFTW3 not found. Please install FFTW3 (e.g., sudo apt install libfftw3-dev)")
endif()

# Find Python and NumPy only if ENABLE_PLOTTING is ON
if(ENABLE_PLOTTING)
    find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
    if(NOT Python3_FOUND)
        message(FATAL_ERROR "Python3/Python3-dev not found. Please install Python3/Python3-dev")
    endif()
endif()

# Print version information
message(STATUS "Found HDF5 version: ${HDF5_VERSION}")
message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "Found FFTW3: ${FFTW_VERSION}")
if(ENABLE_PLOTTING)
    message(STATUS "Plotting enabled")
    message(STATUS "Found Python version: ${Python3_VERSION}")
    message(STATUS "Python executable: ${Python3_EXECUTABLE}")
    message(STATUS "Python include dirs: ${Python3_INCLUDE_DIRS}")
    message(STATUS "NumPy include dirs: ${Python3_NumPy_INCLUDE_DIRS}")
else()
    message(STATUS "Plotting disabled")
endif()

# Include directories
include_directories(
    include
    linearalgebra
    cross_section
    ${HDF5_INCLUDE_DIRS}
)

# Add Python and NumPy include directories only if plotting is enabled
if(ENABLE_PLOTTING)
    include_directories(
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
    )
endif()

# Executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Define ENABLE_PLOTTING macro if plotting is enabled
if(ENABLE_PLOTTING)
    add_definitions(-DENABLE_PLOTTING)
endif()

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    Threads::Threads
    m
    stdc++
    PkgConfig::FFTW
)

# Link Python and NumPy libraries only if plotting is enabled
if(ENABLE_PLOTTING)
    target_link_libraries(${PROJECT_NAME} PRIVATE
        ${Python3_LIBRARIES}
        ${Python3_NumPy_LIBRARIES}
    )
endif()

# Ensure includes are properly set for the target
target_include_directories(${PROJECT_NAME} PRIVATE
    ${HDF5_INCLUDE_DIRS}
)

# Add Python and NumPy include directories for the target only if plotting is enabled
if(ENABLE_PLOTTING)
    target_include_directories(${PROJECT_NAME} PRIVATE
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
    )
endif()

# Test code
enable_testing()

# Test executable for solver performance
add_executable(test_solvers tests/test_solvers.cpp)
target_link_libraries(test_solvers PRIVATE
    Threads::Threads
    m
    PkgConfig::FFTW
)

# Link Python and NumPy for test_solvers only if plotting is enabled
if(ENABLE_PLOTTING)
    target_link_libraries(test_solvers PRIVATE
        ${Python3_LIBRARIES}
        Python3::NumPy
    )
endif()

target_include_directories(test_solvers PRIVATE
    include
    linearalgebra
)

# Add Python and NumPy include directories for test_solvers only if plotting is enabled
if(ENABLE_PLOTTING)
    target_include_directories(test_solvers PRIVATE
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
    )
endif()

# Add test to CTest
add_test(NAME SolverPerformanceTest COMMAND test_solvers)

# Custom target to run tests
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS test_solvers
    COMMENT "Running solver performance tests"
)

# Custom target to clean all build artifacts
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/cmake_install.cmake
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/Makefile
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/${PROJECT_NAME}
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/test_solvers
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaned build directory"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Removing all build artifacts"
)