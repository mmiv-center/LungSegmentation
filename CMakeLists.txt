# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8.9)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

#set(ITK_DIR /Users/hauke/src/itk/InsightToolkit-4.8.1/bin)
set(ITK_DIR /Users/hauke/src/itk/InsightToolkit-5.0.0/bin)
#SET(CMAKE_CXX_FLAGS "-std=c++0x")
SET(CMAKE_CXX_FLAGS "-std=c++11")

project(LungSegmentation)

  FIND_PACKAGE(ITK)
    IF(ITK_FOUND)
        INCLUDE(${ITK_USE_FILE})
    ELSE(ITK_FOUND)
         MESSAGE(FATAL_ERROR "Cannot build ITKApps without ITK.  Please set ITK_DIR.")
    ENDIF(ITK_FOUND)


# Boost
set(boost_min_ver 1.50.0)
set(boost_libs system filesystem timer)
find_package(Boost ${boost_min_ver})

if(Boost_FOUND)
    find_package(Boost ${boost_min_ver} COMPONENTS ${boost_libs})
endif()

#target_link_libraries(your_target ${Boost_LIBRARIES})

add_executable(LungSegmentation LungSegmentation.cxx )
target_link_libraries(LungSegmentation  ${ITK_LIBRARIES} ${Boost_LIBRARIES})

add_executable(MaskRegistration MaskRegistration.cxx )
target_link_libraries(MaskRegistration  ${ITK_LIBRARIES} ${Boost_LIBRARIES})

add_executable(Closing Closing.cxx )
target_link_libraries(Closing  ${ITK_LIBRARIES} ${Boost_LIBRARIES})
