include(ExternalProject)

ExternalProject_Add(geographiclib-external
  GIT_REPOSITORY https://git.code.sf.net/p/geographiclib/code
  GIT_TAG v1.48
  #URL https://downloads.sourceforge.net/project/geographiclib/distrib/GeographicLib-1.48.tar.gz
  CMAKE_ARGS -DGEOGRAPHICLIB_LIB_TYPE=STATIC -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(geographiclib-external BINARY_DIR)
ExternalProject_Get_Property(geographiclib-external SOURCE_DIR)
add_library(geographiclib-lib STATIC IMPORTED)
add_dependencies(geographiclib-lib geographiclib-external)
file(MAKE_DIRECTORY "${SOURCE_DIR}/include")
set_target_properties(geographiclib-lib PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include")
if (WIN32)
set_target_properties(geographiclib-lib PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/lib/Release/${CMAKE_STATIC_LIBRARY_PREFIX}Geographic${CMAKE_STATIC_LIBRARY_SUFFIX}")
else()
set_target_properties(geographiclib-lib PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/src/${CMAKE_STATIC_LIBRARY_PREFIX}Geographic${CMAKE_STATIC_LIBRARY_SUFFIX}")
endif()
