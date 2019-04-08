# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8.9)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

set(ITK_DIR /Users/hauke/src/itk/InsightToolkit-4.8.1/bin)
SET(CMAKE_CXX_FLAGS "-std=c++0x")

project(LungSegmentation)

  FIND_PACKAGE(ITK)
    IF(ITK_FOUND)
        INCLUDE(${ITK_USE_FILE})
    ELSE(ITK_FOUND)
         MESSAGE(FATAL_ERROR "Cannot build ITKApps without ITK.  Please set ITK_DIR.")
    ENDIF(ITK_FOUND)


add_executable(LungSegmentation LungSegmentation.cxx )
target_link_libraries(LungSegmentation  ${ITK_LIBRARIES})

