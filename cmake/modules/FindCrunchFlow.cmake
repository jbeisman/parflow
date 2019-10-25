#.rst:
# FindCrunchFlow
# --------
#
# Find CrunchFlow library
#
# This module finds an installed CrunchFlow library.
#
# This module sets the following variables:
#
# ::
#
#   CRUNCH_FOUND - set to true if a CRUNCH library is found
#   CRUNCH_INCLUDE_DIR - the CRUNCH include directory
#   CRUNCH_LIBRARIES - the CRUNCH libraries
#

include(FindPackageHandleStandardArgs)

if(NOT CRUNCH_ROOT)
    set(CRUNCH_ROOT $ENV{CRUNCH_ROOT})
endif()

set(CRUNCH_INCLUDE_DIR ${CRUNCH_ROOT}/source)

find_library(CRUNCH_LIBRARY NAMES crunchchem
  PATHS ${CRUNCH_INCLUDE_DIR})

set(CRUNCH_INCLUDE_DIRS ${CRUNCH_INCLUDE_DIR})
set(CRUNCH_LIBRARIES ${CRUNCH_LIBRARY})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(CRUNCH DEFAULT_MSG CRUNCH_LIBRARIES CRUNCH_INCLUDE_DIRS)

MARK_AS_ADVANCED(CRUNCH_INCLUDE_DIRS CRUNCH_LIBRARIES)


