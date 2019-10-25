#.rst:
# FindPFLOTRAN
# --------
#
# Find PFLOTRAN library
#
# This module finds an installed PFLOTRAN library.
#
# This module sets the following variables:
#
# ::
#
#   PFLOTRAN_FOUND - set to true if a PFLOTRAN library is found
#   PFLOTRAN_INCLUDE_DIR - the PFLOTRAN include directory
#   PFLOTRAN_LIBRARIES - the PFLOTRAN libraries
#

include(FindPackageHandleStandardArgs)

if(NOT PFLOTRAN_ROOT)
    set(PFLOTRAN_ROOT $ENV{PFLOTRAN_ROOT})
endif()

set(PFLOTRAN_INCLUDE_DIR ${PFLOTRAN_ROOT}/src/pflotran)

find_library(PFLOTRAN_LIBRARY NAMES pflotranchem
  PATHS ${PFLOTRAN_INCLUDE_DIR})

set(PFLOTRAN_INCLUDE_DIRS ${PFLOTRAN_INCLUDE_DIR})
set(PFLOTRAN_LIBRARIES ${PFLOTRAN_LIBRARY})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(PFLOTRAN DEFAULT_MSG PFLOTRAN_LIBRARIES PFLOTRAN_INCLUDE_DIRS)

MARK_AS_ADVANCED(PFLOTRAN_INCLUDE_DIRS PFLOTRAN_LIBRARIES)
