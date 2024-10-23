#---------------------------------------------------
# OpenMM Plumed Plugin
#----------------------------------------------------
SET (CMAKE_C_COMPILER gcc)
SET (CMAKE_CPP_COMPILER g++)
PROJECT(OpenMMPlumed)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Set installation prefix to current directory/install
SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install" CACHE PATH "Installation prefix" FORCE)

# OpenMM directory set relative to current directory
SET(OPENMM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/openmm" CACHE PATH "Where OpenMM is installed")
INCLUDE_DIRECTORIES("${OPENMM_DIR}/include")
LINK_DIRECTORIES("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

# MPI4PY directory with exact path
SET(MPI4PY_DIR "/usr/local/lib/python3.10/site-packages/mpi4py/include" CACHE PATH "Where mpi4py is installed")
INCLUDE_DIRECTORIES("${MPI4PY_DIR}")

# Plumed directories set relative to current directory
SET(PLUMED_LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/plumed/lib" CACHE PATH "Where the Plumed libraries are installed")
SET(PLUMED_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/plumed/include" CACHE PATH "Where the Plumed header files are installed")
INCLUDE_DIRECTORIES("${PLUMED_INCLUDE_DIR}")
LINK_DIRECTORIES("${PLUMED_LIBRARY_DIR}")

# Specify the C++ version we are building for.
SET (CMAKE_CXX_STANDARD 11)

# Set flags for linking on mac
IF(APPLE)
    SET (CMAKE_INSTALL_NAME_DIR "@rpath")
    SET(EXTRA_COMPILE_FLAGS "-msse2 -stdlib=libc++")
ENDIF(APPLE)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(PLUMED_PLUGIN_SOURCE_SUBDIRS openmmapi)

# Set the library name
SET(PLUMED_LIBRARY_NAME OpenMMPlumed)
SET(SHARED_PLUMED_TARGET ${PLUMED_LIBRARY_NAME})

# These are all the places to search for header files which are to be part of the API.
SET(API_INCLUDE_DIRS "openmmapi/include" "openmmapi/include/internal")

# Locate header files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
ENDFOREACH(dir)

# Collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)
FOREACH(subdir ${PLUMED_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# Build MPI and make it required
FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${MPI_INCLUDE_PATH})

# Create the library.
ADD_LIBRARY(${SHARED_PLUMED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_PLUMED_TARGET}
    PROPERTIES COMPILE_FLAGS "-DPLUMED_BUILDING_SHARED_LIBRARY ${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_COMPILE_FLAGS}")
TARGET_LINK_LIBRARIES(${SHARED_PLUMED_TARGET} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} OpenMM plumed)
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_PLUMED_TARGET})

# install headers
FILE(GLOB API_ONLY_INCLUDE_FILES "openmmapi/include/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES} DESTINATION include)
FILE(GLOB API_ONLY_INCLUDE_FILES_INTERNAL "openmmapi/include/internal/*.h")
INSTALL (FILES ${API_ONLY_INCLUDE_FILES_INTERNAL} DESTINATION include/internal)

# Build the implementations for different platforms
ADD_SUBDIRECTORY(platforms/reference)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
FIND_PACKAGE(OpenCL QUIET)
IF(OPENCL_FOUND)
    SET(PLUMED_BUILD_OPENCL_LIB OFF CACHE BOOL "Build implementation for OpenCL")
ELSE(OPENCL_FOUND)
    SET(PLUMED_BUILD_OPENCL_LIB OFF CACHE BOOL "Build implementation for OpenCL")
ENDIF(OPENCL_FOUND)
IF(PLUMED_BUILD_OPENCL_LIB)
    ADD_SUBDIRECTORY(platforms/opencl)
ENDIF(PLUMED_BUILD_OPENCL_LIB)

FIND_PACKAGE(CUDA QUIET)
IF(CUDA_FOUND)
    SET(PLUMED_BUILD_CUDA_LIB OFF CACHE BOOL "Build implementation for CUDA")
ELSE(CUDA_FOUND)
    SET(PLUMED_BUILD_CUDA_LIB OFF CACHE BOOL "Build implementation for CUDA")
ENDIF(CUDA_FOUND)
IF(PLUMED_BUILD_CUDA_LIB)
    ADD_SUBDIRECTORY(platforms/cuda)
ENDIF(PLUMED_BUILD_CUDA_LIB)

# Build the Python API
FIND_PROGRAM(PYTHON_EXECUTABLE python)
FIND_PROGRAM(SWIG_EXECUTABLE swig)
IF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(PLUMED_BUILD_PYTHON_WRAPPERS ON CACHE BOOL "Build wrappers for Python")
ELSE(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
    SET(PLUMED_BUILD_PYTHON_WRAPPERS OFF CACHE BOOL "Build wrappers for Python")
ENDIF(PYTHON_EXECUTABLE AND SWIG_EXECUTABLE)
IF(PLUMED_BUILD_PYTHON_WRAPPERS)
    ADD_SUBDIRECTORY(python)
ENDIF(PLUMED_BUILD_PYTHON_WRAPPERS)