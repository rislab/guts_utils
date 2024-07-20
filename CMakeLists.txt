cmake_minimum_required(VERSION 2.8.3)

project(guts_utils)
set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)

include_directories(include
  ${EIGEN3_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
  SHARED
  src/GUTSUtils.cc
)

target_link_libraries(${PROJECT_NAME}
  PRIVATE Eigen3::Eigen
)

set_target_properties(${PROJECT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX})

install(TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION "include/${PROJECT_NAME}"
  FILES_MATCHING PATTERN "*.h"
)

install(DIRECTORY cmake/
  DESTINATION "share/${PROJECT_NAME}"
  FILES_MATCHING PATTERN "*.cmake"
)
