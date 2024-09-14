# Find the GMP header and library
find_path(GMP_INCLUDES
  NAMES gmp.h
  PATHS /usr/include /usr/local/include
)

find_library(GMP_LIBRARIES
  NAMES gmp
  PATHS /usr/lib /usr/local/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDES GMP_LIBRARIES)

# Define the imported target for GMP
if (GMP_FOUND)
  add_library(GMP::GMP INTERFACE IMPORTED)
  set_target_properties(GMP::GMP PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDES}"
    INTERFACE_LINK_LIBRARIES "${GMP_LIBRARIES}"
  )
endif()
