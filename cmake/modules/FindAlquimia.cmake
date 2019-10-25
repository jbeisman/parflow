#.rst:
# FindAlquimia
# --------
#
# Find Alquimia library
#
# This module finds an installed Alquimia library.
#
# This module sets the following variables:
#
# ::
#
#   ALQUIMIA_FOUND - set to true if a ALQUIMIA library is found
#   ALQUIMIA_INCLUDE_DIR - the ALQUIMIA include directory
#   ALQUIMIA_LIBRARIES - the ALQUIMIA libraries
#

include(FindPackageHandleStandardArgs)

if(NOT ALQUIMIA_ROOT)
    set(ALQUIMIA_ROOT $ENV{ALQUIMIA_ROOT})
endif()

set(ALQUIMIA_INCLUDE_DIR ${ALQUIMIA_ROOT}/include)

find_library(ALQUIMIA_LIBRARY NAMES alquimia
  PATHS ${ALQUIMIA_ROOT}/lib)

set(ALQUIMIA_INCLUDE_DIRS ${ALQUIMIA_INCLUDE_DIR})
set(ALQUIMIA_LIBRARIES ${ALQUIMIA_LIBRARY})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(ALQUIMIA DEFAULT_MSG ALQUIMIA_LIBRARIES ALQUIMIA_INCLUDE_DIRS)

MARK_AS_ADVANCED(ALQUIMIA_INCLUDE_DIRS ALQUIMIA_LIBRARIES)


