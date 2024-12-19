
# CMakeLists.txt

cmake_minimum_required(VERSION 3.15)
project(XSol-TraderAI VERSION 1.0.0 LANGUAGES CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

# Find required packages
find_package(CUDA REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem)
find_package(OpenSSL REQUIRED)
find_package(CURL REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${CUDA_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${OPENSSL_INCLUDE_DIR}
    ${CURL_INCLUDE_DIRS}
)

# Source files
file(GLOB_RECURSE SOURCES
    "src/*.cpp"
    "src/*.cu"
    "src/models/*.cpp"
    "src/utils/*.cpp"
)

# Create main executable
add_executable(xsol_trader_ai ${SOURCES})

# Link libraries
target_link_libraries(xsol_trader_ai PRIVATE
    ${CUDA_LIBRARIES}
    ${Boost_LIBRARIES}
    ${OPENSSL_LIBRARIES}
    ${CURL_LIBRARIES}
    nlohmann_json::nlohmann_json
    pthread
)

# Enable CUDA support
set_target_properties(xsol_trader_ai PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Install targets
install(TARGETS xsol_trader_ai
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Testing
enable_testing()
add_subdirectory(tests)

# Documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_subdirectory(docs)
endif()

# Package configuration
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/XSolTraderAIConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Install configuration
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/XSolTraderAIConfigVersion.cmake"
    DESTINATION lib/cmake/XSolTraderAI
)

# Custom targets
add_custom_target(format
    COMMAND clang-format -i ${SOURCES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)