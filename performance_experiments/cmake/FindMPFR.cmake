# Find the MPFR header and library
find_path(MPFR_INCLUDES
  NAMES mpfr.h
  PATHS /usr/include /usr/local/include
)

find_library(MPFR_LIBRARIES
  NAMES mpfr
  PATHS /usr/lib /usr/local/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_INCLUDES MPFR_LIBRARIES)

# Define the imported target for MPFR
if (MPFR_FOUND)
  add_library(MPFR::MPFR INTERFACE IMPORTED)
  set_target_properties(MPFR::MPFR PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDES}"
    INTERFACE_LINK_LIBRARIES "${MPFR_LIBRARIES}"
  )
endif()
