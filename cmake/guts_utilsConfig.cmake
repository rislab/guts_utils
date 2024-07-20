# Finds guts_utils header files

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)
set(guts_utils_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/guts_utils/include")
if (APPLE)
  set(guts_utils_LIBRARIES "${PACKAGE_PREFIX_DIR}/guts_utils/libguts_utils.dylib")
elseif(UNIX)
  set(guts_utils_LIBRARIES "${PACKAGE_PREFIX_DIR}/guts_utils/libguts_utils.so")
else()
  message(FATAL_ERROR "Operating system not supported")
endif()
set(guts_utils_LIBRARY_DIR "${PACKAGE_PREFIX_DIR}/guts_utils")
