project (grav)
cmake_minimum_required (VERSION 3.1.0)

find_package (GeographicLib REQUIRED)

add_executable (${PROJECT_NAME} NormalGravity.cpp)
target_link_libraries (${PROJECT_NAME} ${GeographicLib_LIBRARIES})
