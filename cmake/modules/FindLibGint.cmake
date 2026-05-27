#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2025- ETH Zurich
#
# author: Rocco Meli

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_find_libraries(LIBGINT ("Gint")
cp2k_include_dirs(LIBGINT "libgint_wrapper.mod")

if(CP2K_LIBGINT_INCLUDE_DIRS)
  find_package_handle_standard_args(
    libgint DEFAULT_MSG CP2K_LIBGINT_LINK_LIBRARIES
    CP2K_LIBGINT_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(libgint DEFAULT_MSG
                                    CP2K_LIBGINT_LINK_LIBRARIES)
endif()

if(NOT TARGET cp2k::libGint::libGint)
  add_library(cp2k::libGint::libGint INTERFACE IMPORTED)
  set_target_properties(
    cp2k::libGint::libGint PROPERTIES INTERFACE_LINK_LIBRARIES
                                         "${CP2K_LIBSGINT_LINK_LIBRARIES}")
  if(CP2K_LIBGINT_INCLUDE_DIRS)
    set_target_properties(
      cp2k::libGint::libGint PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                           "${CP2K_LIBGINT_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_LIBGINT_ROOT CP2K_LIBGINT_INCLUDE_DIRS
                 CP2K_LIBGINT_LINK_LIBRARIES CP2K_LIBGINT_LIBRARIES)
