#----------------------------------------------------------------------------
# Setup the project
#
project(FieldCalibration C CXX)

#----------------------------------------------------------------------------
# cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
#
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

#----------------------------------------------------------------------------
# Use C++11 extention and O3 optimization (-founding-math because of CGAL)
#
set(CMAKE_CXX_FLAGS "-O3 -frounding-math")
add_definitions("-std=c++11")

#----------------------------------------------------------------------------
# Find CGAL library
#
find_package(CGAL QUIET)
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED)
add_definitions(${ROOT_CXX_FLAGS})

#----------------------------------------------------------------------------
# Find package for multi-threading
#
find_package (Threads)

#----------------------------------------------------------------------------
# Add sources and headers
#
set(SOURCE_FILES source/Laser.cpp source/LaserTrack.cpp source/Laser.cpp source/Matrix3x3.cpp source/TPCVolumeHandler source/Interpolation3D.cpp)
set(HEADER_FILES include/Laser.hpp include/LaserTrack.hpp include/Laser.hpp include/ThreeVector.hpp include/Matrix3x3.hpp include/TPCVolumeHandler include/Interpolation3D.hpp)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${ROOT_INCLUDE_DIR})


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(FieldCal FieldCal.cpp ${HEADER_FILES} ${SOURCE_FILES} )
target_link_libraries(FieldCal ${ROOT_LIBRARIES} ${CGAL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

#----------------------------------------------------------------------------
# Copy all required files to build directory.
#
#set(FieldCal_Files Field.root)

#foreach(_File ${FieldCal_Files})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/${_File}
#    ${PROJECT_BINARY_DIR}/${_File}
#    COPYONLY
#    )
#endforeach()


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS FieldCal DESTINATION bin)

