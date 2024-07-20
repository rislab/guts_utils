# Finds guts_utils header files

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)
set(guts_utils_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/guts_utils/include")
set(guts_utils_LIBRARIES "${PACKAGE_PREFIX_DIR}/guts_utils/libguts_utils.dylib")
set(guts_utils_LIBRARY_DIR "${PACKAGE_PREFIX_DIR}/guts_utils")
