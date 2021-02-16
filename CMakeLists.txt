cmake_minimum_required(VERSION 3.1)
project(geographiclib_conversions)

# Add module path for geographiclib
find_package(catkin QUIET)
if (catkin_FOUND)
  set(USE_ROS TRUE)
endif()

include(3rdparty/External_geographiclib.cmake)

if (USE_ROS)
  find_package(catkin REQUIRED COMPONENTS roscpp geometry_msgs std_msgs)
  catkin_package(
    INCLUDE_DIRS include 
    LIBRARIES ${PROJECT_NAME}
    CATKIN_DEPENDS 
  )
endif()

## catkin specific configuration

add_library(${PROJECT_NAME} 
  src/geographiclib_conversions.cpp
)
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} geographiclib-lib)
target_include_directories(${PROJECT_NAME} PUBLIC include ${catkin_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}_main test/geographiclib_main.cpp) 
target_link_libraries(${PROJECT_NAME}_main PUBLIC ${PROJECT_NAME})

if (USE_ROS)

find_package(rostest REQUIRED)

if(CATKIN_ENABLE_TESTING)
  find_package(rostest REQUIRED)
  catkin_add_gtest(${PROJECT_NAME}_test test/geographiclib_test.cpp)
  target_link_libraries(${PROJECT_NAME}_test ${catkin_LIBRARIES} ${PROJECT_NAME})
endif()

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})

install(TARGETS ${PROJECT_NAME}
  DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
endif()
