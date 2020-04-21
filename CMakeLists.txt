cmake_minimum_required(VERSION 3.13)

project(EnDAS VERSION 0.1)

include (GNUInstallDirs)
include (WriteCompilerDetectionHeader)


if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed. Create 'build' sub-directory and run CMake from there: \n> mkdir build\n> cd build\n> cmake ..")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(DISABLE_BLAS_LAPACK "Disable usage of BLAS/LAPACK (Eigen's internal implementation will be used)" OFF)
option(DISABLE_PLOTTING "Disable plotting in examples" OFF)
option(DISABLE_PROFILING "Disable collection of performance profiling information in algorithms" OFF)


find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if (NOT DISABLE_BLAS_LAPACK)
   find_package(BLAS)    
   find_package(LAPACK)    
   find_package(LAPACKE)    
endif()

if (NOT DISABLE_PLOTTING)
   find_package(Python2 COMPONENTS Development)
endif()


set (ENDAS_USING_BLAS 0)
set (ENDAS_USING_LAPACK 0)
set (ENDAS_USING_PYTHON 0)
set (ENDAS_PROFILING_DISABLED 0)


message(STATUS "Configuration summary:")
message(STATUS "-------------------------------------------------------")

message(STATUS "Prefix   : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Eigen3   : ${EIGEN3_INCLUDE_DIR}")

if (BLAS_FOUND)
   message(STATUS "BLAS     : yes (${BLAS_LIBRARIES})")
   set (ENDAS_USING_BLAS 1)
else()
   message(STATUS "BLAS     : no")
endif()

if (LAPACK_FOUND AND LAPACKE_FOUND)
   message(STATUS "LAPACK   : yes (${LAPACKE_LIBRARIES})")
   set (ENDAS_USING_LAPACK 1)
else()
   message(STATUS "LAPACK   : no")
endif()

if (TARGET Python2::Python)
   get_target_property(PYTHON2_VERSION Python2::Python VERSION)
   message(STATUS "Python2  : yes (${Python2_LIBRARIES})")
   set (ENDAS_USING_PYTHON 1)
else()
   message(STATUS "Python2  : no")
endif()

if (DISABLE_PROFILING)
   message(STATUS "Profiling: disabled")
   set (ENDAS_PROFILING_DISABLED 1)
else()
   message(STATUS "Profiling: enabled")
endif()

message(STATUS "-------------------------------------------------------")


write_compiler_detection_header(
   FILE ${PROJECT_SOURCE_DIR}/include/Endas/Compiler.h
   PREFIX ENDAS
   COMPILERS GNU Clang MSVC Intel
   FEATURES cxx_variadic_templates 
)

set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/thirdparty)


configure_file(${PROJECT_SOURCE_DIR}/include/Endas/Config.h.in ${PROJECT_SOURCE_DIR}/include/Endas/Config.h)


add_subdirectory(src)
add_subdirectory(examples)




