cmake_minimum_required(VERSION 3.22)

# Read version from root VERSION file (Task A20)
file(READ "${CMAKE_SOURCE_DIR}/VERSION" ELLE_VERSION_STR)
string(STRIP "${ELLE_VERSION_STR}" ELLE_VERSION_STR)
project(ELLE_GeoSim VERSION "${ELLE_VERSION_STR}" LANGUAGES C CXX)
add_compile_definitions(ELLE_VERSION="${PROJECT_VERSION}")
# Legacy code relies on deprecated Boost progress/timer headers. Define the
# opt‑in macro globally so individual targets don’t fail on inclusion.
add_compile_definitions(BOOST_TIMER_ENABLE_DEPRECATED=1)

# Conan (opt-in) acknowledgement
if(DEFINED CMAKE_PROJECT_TOP_LEVEL_INCLUDES)
  message(STATUS "Conan toolchain active")
endif()

# Options
option(ELLE_WARN_AS_ERROR "Treat warnings as errors" OFF)
option(ELLE_HEADLESS "Build without GUI/display backends" ON)

# New optional features (Task A12)
option(ELLE_BUILD_GUI "Build GUI components with wxWidgets/X11" OFF)
option(ELLE_USE_MPI "Enable MPI features" OFF)
option(ELLE_BUILD_ALL_PROCESSES "Build all legacy process executables" OFF)
set(ELLE_PROCESS "" CACHE STRING "Build only the named legacy process (stem). Ignored if ELLE_BUILD_ALL_PROCESSES=ON")

# Language standards and compile database
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Install dir helpers and basic vendor/contact metadata (Task A19)
include(GNUInstallDirs)
set(ELLE_VENDOR "ELLE Project")
set(CPACK_PACKAGE_CONTACT "maintainers@example.org")

# Build type default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Warnings
if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
  if(ELLE_WARN_AS_ERROR)
    add_compile_options(-Werror)
  endif()
endif()

# Common dependencies (enable as Analysis decides)
# find_package(Eigen3 3.3 QUIET)
# find_package(HDF5 COMPONENTS C HL QUIET)
# find_package(BLAS QUIET)
# find_package(LAPACK QUIET)

# MPI (optional, guarded)
if(ELLE_USE_MPI)
  find_package(MPI QUIET COMPONENTS C CXX)
endif()

# ------------------------------------------------------------
# Headless core dependencies (no GUI/wx/X11 for first pass)
# ------------------------------------------------------------
find_package(ZLIB REQUIRED)
# GSL provides both GSL::gsl and GSL::gslcblas
find_package(GSL REQUIRED)
# IPhreeqc (via Conan if available). Optional: if not found, the code falls back to local stub.
find_package(IPhreeqc QUIET)
# Boost (required for string algorithms in latte/phreeqc)
find_package(Boost REQUIRED)
# OpenMP (needed by phreeqc-based processes for omp_get_max_threads). Optional.
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found: enabling parallel runtime for phreeqc processes")
else()
  message(STATUS "OpenMP not found: phreeqc processes will run single-threaded; consider installing libgomp / compiler OpenMP support")
endif()

# ------------------------------------------------------------
# Legacy core sources (elle/elle/basecode)
# ------------------------------------------------------------
set(ELLE_BASECODE_DIR "${CMAKE_SOURCE_DIR}/elle/elle/basecode")
set(ELLE_STATSCODE_DIR "${CMAKE_SOURCE_DIR}/elle/elle/statscode") # contains stats.h
set(ELLE_GPC_DIR       "${CMAKE_SOURCE_DIR}/elle/elle/utilities/gpc")
set(ELLE_MATRIX_DIR    "${CMAKE_SOURCE_DIR}/elle/elle/utilities/matrix")

# Collect C++ and C sources from the legacy basecode directory
file(GLOB_RECURSE ELLE_CORE_CXX CONFIGURE_DEPENDS
     "${ELLE_BASECODE_DIR}/*.cc" "${ELLE_BASECODE_DIR}/*.cpp")
file(GLOB_RECURSE ELLE_CORE_C   CONFIGURE_DEPENDS
     "${ELLE_BASECODE_DIR}/*.c")

# Exclude obvious GUI/wx/X-specific files if any slipped in
list(FILTER ELLE_CORE_CXX EXCLUDE REGEX ".*/wx.*\\.(cc|cpp)$")
list(FILTER ELLE_CORE_CXX EXCLUDE REGEX ".*/display_x\\.(cc|cpp)$")
# Keep GUI/wx out, but include core math/utils needed by basecode

set(ELLE_CORE_UTIL_CXX
  "${ELLE_MATRIX_DIR}/mat.cc"         # matrix/euler helpers
  "${ELLE_GPC_DIR}/gpcclip.cc"       # polygon clipping wrapper
)
set(ELLE_CORE_UTIL_C
  "${ELLE_GPC_DIR}/gpc.c"            # GPC implementation (C)
)

add_library(elle_core STATIC ${ELLE_CORE_CXX} ${ELLE_CORE_C} ${ELLE_CORE_UTIL_CXX} ${ELLE_CORE_UTIL_C} "${ELLE_STATSCODE_DIR}/statsfns.cc")

# Ensure Triangle builds as a library (provides triangulate())
set(ELLE_TRIANGLE_C "${ELLE_BASECODE_DIR}/triangle.c")
if(EXISTS "${ELLE_TRIANGLE_C}")
  set_source_files_properties("${ELLE_TRIANGLE_C}" PROPERTIES COMPILE_DEFINITIONS TRILIBRARY)
endif()

# Expose exact on-disk include paths required by basecode and its direct deps.
# Specifically, stats.h lives in ${ELLE_STATSCODE_DIR} and must be visible to all TUs.
target_include_directories(elle_core PUBLIC
  "${ELLE_BASECODE_DIR}"
  "${ELLE_MATRIX_DIR}"
  "${ELLE_STATSCODE_DIR}" # for stats.h and related stats APIs
  "${ELLE_GPC_DIR}"       # for gpcclip.h (polygon clipping)
  "${CMAKE_SOURCE_DIR}/include" # for common portability macros
)

# Core library is always GUI-free: define NODISP for elle_core
# (GUI functionality is supplied by the elle_gui target instead.)
# modernize: Do not propagate NODISP to consumers (e.g., GUI); keep it internal (analysis fast-track)
target_compile_definitions(elle_core PRIVATE NODISP)


# Note: runopts.cc no longer defines stub Run_App/SetupApp even with NODISP,
# so we keep normal NODISP build for elle_core and avoid duplicate symbols.

# Link core dependencies
target_link_libraries(elle_core
  PUBLIC
    ZLIB::ZLIB
    GSL::gsl
    GSL::gslcblas
)

# If MPI is discovered and later needed, we can link it here conditionally
# if(MPI_CXX_FOUND)
#   target_link_libraries(elle_core PUBLIC MPI::MPI_CXX)
#   target_compile_definitions(elle_core PUBLIC ELLE_USE_MPI=1)
# endif()

# Optional LAPACK (Task A12)
find_package(LAPACK QUIET)
if(LAPACK_FOUND)
  target_link_libraries(elle_core PUBLIC LAPACK::LAPACK)
endif()

# If MPI requested and found, link and define
if(ELLE_USE_MPI)
  if(MPI_FOUND)
    target_compile_definitions(elle_core PUBLIC ELLE_USE_MPI)
    # Link both C and C++ MPI components when available
    if(TARGET MPI::MPI_CXX)
      target_link_libraries(elle_core PUBLIC MPI::MPI_CXX)
    endif()
    if(TARGET MPI::MPI_C)
      target_link_libraries(elle_core PUBLIC MPI::MPI_C)
    endif()
    # Optionally add MPI sources if present (mpi_*.cc)
    file(GLOB_RECURSE ELLE_MPI_SOURCES CONFIGURE_DEPENDS
         "${ELLE_BASECODE_DIR}/mpi_*.cc" "${CMAKE_SOURCE_DIR}/elle/**/mpi_*.cc")
    if(ELLE_MPI_SOURCES)
      target_sources(elle_core PRIVATE ${ELLE_MPI_SOURCES})
    endif()
  else()
    message(STATUS "MPI disabled: not found")
  endif()
endif()

# ------------------------------------------------------------
# Minimal headless CLI for link/smoke tests (behavior-neutral)
# ------------------------------------------------------------
add_executable(elle_cli src/elle_cli.cc)
target_link_libraries(elle_cli PRIVATE elle_core)

if(UNIX)
  # Link libm if needed for transitive math symbols
  target_link_libraries(elle_cli PRIVATE m)
endif()

# Install rules (Task A19)
install(TARGETS elle_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# ------------------------------------------------------------
# Optional GUI target (wxWidgets) — guarded (Task A12)
# ------------------------------------------------------------
if(ELLE_BUILD_GUI)
  # modernize: Prefer minimal wxWidgets components for initial triage (analysis fast-track)
  # Need html component for wxHtmlWindow used by showelle
  find_package(wxWidgets 3.2 QUIET COMPONENTS core base html)
  if(wxWidgets_FOUND)
    # Build a GUI-enabled core variant without NODISP so StartApp can call Run_App
    add_library(elle_core_gui STATIC ${ELLE_CORE_CXX} ${ELLE_CORE_C} ${ELLE_CORE_UTIL_CXX} ${ELLE_CORE_UTIL_C} "${ELLE_STATSCODE_DIR}/statsfns.cc")
    target_include_directories(elle_core_gui PUBLIC
      "${ELLE_BASECODE_DIR}"
      "${ELLE_MATRIX_DIR}"
      "${ELLE_STATSCODE_DIR}"
      "${ELLE_GPC_DIR}"
      "${CMAKE_SOURCE_DIR}/elle/elle/wxplotcode" # for settings.h and GUI headers used by display.cc
      "${CMAKE_SOURCE_DIR}/include")
    target_link_libraries(elle_core_gui
      PUBLIC ZLIB::ZLIB GSL::gsl GSL::gslcblas)
    if(LAPACK_FOUND)
      target_link_libraries(elle_core_gui PUBLIC LAPACK::LAPACK)
    endif()
    if(ELLE_USE_MPI)
      if(MPI_FOUND)
        target_compile_definitions(elle_core_gui PUBLIC ELLE_USE_MPI)
        if(TARGET MPI::MPI_CXX)
          target_link_libraries(elle_core_gui PUBLIC MPI::MPI_CXX)
        endif()
        if(TARGET MPI::MPI_C)
          target_link_libraries(elle_core_gui PUBLIC MPI::MPI_C)
        endif()
        file(GLOB_RECURSE ELLE_MPI_SOURCES CONFIGURE_DEPENDS
             "${ELLE_BASECODE_DIR}/mpi_*.cc" "${CMAKE_SOURCE_DIR}/elle/**/mpi_*.cc")
        if(ELLE_MPI_SOURCES)
          target_sources(elle_core_gui PRIVATE ${ELLE_MPI_SOURCES})
        endif()
      endif()
    endif()
    # Collect potential GUI sources from legacy tree
  file(GLOB_RECURSE ELLE_GUI_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_SOURCE_DIR}/elle/elle/wxplotcode/*.cc"  # all GUI sources
     "${CMAKE_SOURCE_DIR}/elle/**/wx*.cc"            # any wx*-prefixed files elsewhere
     "${CMAKE_SOURCE_DIR}/elle/**/display_x*.cc")
    if(ELLE_GUI_SOURCES)
  add_executable(elle_gui ${ELLE_GUI_SOURCES})
  # Link GUI app against GUI-enabled core so display.cc and friends are included
  target_link_libraries(elle_gui PRIVATE elle_core_gui wxWidgets::wxWidgets)
      if(UNIX)
        target_link_libraries(elle_gui PRIVATE m)
      endif()
      # Install only if target exists
      install(TARGETS elle_gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  # latte-clean-old with GUI: reuse latte sources and link GUI + wx (headless stubs removed)
  set(LATTE_OLD_GUI_SOURCES ${LATTE_OLD_SOURCES})
  add_executable(latte_clean_old_gui ${LATTE_OLD_GUI_SOURCES} ${ELLE_GUI_SOURCES})
      target_include_directories(latte_clean_old_gui PRIVATE
        "${LATTE_OLD_DIR}"
        "${ELLE_BASECODE_DIR}"
        "${ELLE_MATRIX_DIR}"
        "${ELLE_STATSCODE_DIR}"
        "${ELLE_GPC_DIR}"
        "${CMAKE_SOURCE_DIR}/include")
      target_link_libraries(latte_clean_old_gui PRIVATE elle_core_gui wxWidgets::wxWidgets)
      if(UNIX)
        target_link_libraries(latte_clean_old_gui PRIVATE m)
      endif()
      install(TARGETS latte_clean_old_gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    else()
      message(STATUS "GUI disabled: no GUI sources found")
    endif()
  else()
    message(STATUS "GUI disabled: wxWidgets not found")
  endif()
endif()

# ------------------------------------------------------------
# Library install (static) — optional for consumers/tools
# ------------------------------------------------------------
install(TARGETS elle_core
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/elle)

# Note: Public headers are not yet curated. This project is primarily an
# application; skipping header install avoids exposing unstable internals.
# If needed later, group headers under include/elle and install that dir.
#
# install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/elle/
#         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/elle
#         FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

# ------------------------------------------------------------
# RPATH/RUNPATH sanity for install tree (Task A19)
# ------------------------------------------------------------
set(CMAKE_SKIP_BUILD_RPATH OFF)
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set(CMAKE_INSTALL_RPATH "")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)

# ------------------------------------------------------------
# CPack configuration (TGZ + DEB) (Task A19)
# ------------------------------------------------------------
set(CPACK_PACKAGE_NAME "elle")
set(CPACK_PACKAGE_VENDOR "${ELLE_VENDOR}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ELLE Geological Simulation (headless CLI)")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_GENERATOR "TGZ;DEB")

# Debian-specific metadata (safe defaults)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "${CPACK_PACKAGE_CONTACT}")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
# Keep runtime deps light for headless CLI
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libgsl0ldbl | libgsl27, zlib1g")

include(CPack)

# ------------------------------------------------------------
# Process: latte-clean-old (headless build) 
# ------------------------------------------------------------
# Build the legacy latte-clean-old process against elle_core without GUI.
# We provide tiny stubs for SetupApp/Run_App so the old mains link in NODISP.
set(LATTE_OLD_DIR "${CMAKE_SOURCE_DIR}/elle/elle/processes/latte-clean-old")

set(LATTE_OLD_SOURCES
  "${LATTE_OLD_DIR}/elle.main.cc"
  "${LATTE_OLD_DIR}/experiment.cc"
  "${LATTE_OLD_DIR}/lattice.cc"
  "${LATTE_OLD_DIR}/particle.cc"
  "${LATTE_OLD_DIR}/phase_base.cc"
  "${LATTE_OLD_DIR}/phase_lattice.cc"
  "${LATTE_OLD_DIR}/fluid_lattice.cc"
  "${LATTE_OLD_DIR}/heat_lattice.cc"
  "${LATTE_OLD_DIR}/pressure_lattice.cc"
  "${LATTE_OLD_DIR}/graingrowth.cc"
  "${LATTE_OLD_DIR}/min_trans.cc"
  "${LATTE_OLD_DIR}/veingrowth.cc"
)

# Exclude GUI helper and alternate mains not needed for headless build
add_executable(latte_clean_old ${LATTE_OLD_SOURCES})
target_include_directories(latte_clean_old PRIVATE
  "${LATTE_OLD_DIR}"
  "${ELLE_BASECODE_DIR}"
  "${ELLE_MATRIX_DIR}"
  "${ELLE_STATSCODE_DIR}"
  "${ELLE_GPC_DIR}"
  "${CMAKE_SOURCE_DIR}/include"
)
target_link_libraries(latte_clean_old PRIVATE elle_core)
if(UNIX)
  target_link_libraries(latte_clean_old PRIVATE m)
endif()

# ------------------------------------------------------------
# Auto-discovered legacy process executables (selective)
# ------------------------------------------------------------
# Options:
#   -DELLE_BUILD_ALL_PROCESSES=ON  -> build everything discoverable
#   -DELLE_PROCESS=name            -> build only that process (stem) if found (mutually exclusive with ALL)
# If neither set, only explicitly defined processes (latte, latte-clean-old, latte) are built.

if(ELLE_BUILD_ALL_PROCESSES AND ELLE_PROCESS)
  message(WARNING "Both ELLE_BUILD_ALL_PROCESSES and ELLE_PROCESS specified; ELLE_PROCESS will be ignored.")
endif()

if(ELLE_BUILD_ALL_PROCESSES OR ELLE_PROCESS)
  set(ELLE_PROCESSES_ROOT "${CMAKE_SOURCE_DIR}/elle/elle/processes")
  file(GLOB_RECURSE ELLE_PROCESS_MAIN_FILES CONFIGURE_DEPENDS
    "${ELLE_PROCESSES_ROOT}/*.main.cc")
  set(EXCLUDE_PROCESS_DIR_REGEX ".*/elle/elle/processes/(latte|latte-clean-old|latte-noP)/.*")
  set(ADDED_PROCESS_TARGETS "")
  set(NEEDS_FORTRAN FALSE)

  foreach(MAIN_FILE IN LISTS ELLE_PROCESS_MAIN_FILES)
    if(MAIN_FILE MATCHES "${EXCLUDE_PROCESS_DIR_REGEX}")
      continue()
    endif()
    get_filename_component(MAIN_DIR "${MAIN_FILE}" DIRECTORY)
    get_filename_component(MAIN_NAME "${MAIN_FILE}" NAME_WE)
    string(REGEX REPLACE "\\.main$" "" STEM "${MAIN_NAME}")
    if(STEM STREQUAL "main")
      get_filename_component(DIR_STEM "${MAIN_DIR}" NAME)
      set(TARGET_NAME "${DIR_STEM}")
      file(GLOB DIR_CC CONFIGURE_DEPENDS "${MAIN_DIR}/*.cc")
      set(SRCS ${DIR_CC})
    else()
      set(TARGET_NAME "${STEM}")
      file(GLOB STEM_CC CONFIGURE_DEPENDS "${MAIN_DIR}/${STEM}*.cc")
      file(GLOB STEM_C  CONFIGURE_DEPENDS "${MAIN_DIR}/${STEM}*.c")
      file(GLOB STEM_F  CONFIGURE_DEPENDS "${MAIN_DIR}/${STEM}*.f")
      set(SRCS ${STEM_CC} ${STEM_C} ${STEM_F} ${MAIN_FILE})
    endif()
    list(REMOVE_DUPLICATES SRCS)
    if(ELLE_PROCESS AND NOT ELLE_BUILD_ALL_PROCESSES)
      string(TOLOWER "${ELLE_PROCESS}" REQ_LOWER)
      string(TOLOWER "${TARGET_NAME}" TGT_LOWER)
      if(NOT REQ_LOWER STREQUAL TGT_LOWER)
        continue()
      endif()
    endif()
    # Skip if already defined
    if(TARGET ${TARGET_NAME})
      continue()
    endif()
    foreach(_s IN LISTS SRCS)
      if(_s MATCHES "\\.f$")
        # Enable Fortran as soon as a .f source is detected so targets using
        # Fortran (e.g. tbh) pick up the correct linker/runtime flags.
        if(NOT CMAKE_Fortran_COMPILER)
          enable_language(Fortran)
        endif()
        set(NEEDS_FORTRAN TRUE)
      endif()
    endforeach()
    # Process-specific source pruning for legacy overlaps
    if(TARGET_NAME STREQUAL "gbe")
      # Both gbe.elle.cc and gbe_growth_euler.cc implement the same functions.
      # Keep the shorter gbe.elle.cc variant and drop gbe_growth_euler.cc to
      # avoid multiple definition link errors.
      list(FILTER SRCS EXCLUDE REGEX ".*/gbe_growth_euler.cc$")
    endif()
    # Special case: NewFlow_phreeqc-2 process (target name 'elle') needs all
    # directory sources, not just those matching the stem prefix, because its
    # main file is elle.main.cc but supporting implementation files have
    # unrelated basenames (experiment.cc, fluid_lattice.cc, etc.). Also allow
    # building without system IPhreeqc by using the bundled stub header.
    if(TARGET_NAME STREQUAL "elle" AND MAIN_DIR MATCHES ".*/NewFlow_phreeqc-2$")
      file(GLOB _dir_all_cc CONFIGURE_DEPENDS "${MAIN_DIR}/*.cc")
      # Ensure main file present
      list(APPEND _dir_all_cc "${MAIN_DIR}/elle.main.cc")
      # Drop GUI-only source requiring wxWidgets for headless build
      list(FILTER _dir_all_cc EXCLUDE REGEX ".*/startfromgui.cc$")
      set(SRCS ${_dir_all_cc})
      list(REMOVE_DUPLICATES SRCS)
    endif()
    if(TARGET_NAME STREQUAL "jkbmelt")
      # jkbmelt needs helper sources that don't share the jkbmelt* stem
      # (e.g. meltmovenode.cc). Add any melt*.cc present in the same directory.
      get_filename_component(_proc_dir "${MAIN_FILE}" DIRECTORY)
      file(GLOB _extra_melt CONFIGURE_DEPENDS "${_proc_dir}/melt*.cc")
      list(APPEND SRCS ${_extra_melt})
      list(REMOVE_DUPLICATES SRCS)
    endif()
    if(TARGET_NAME STREQUAL "gbm")
      # gbm needs movenode.cc which doesn't share the gbm* stem
      get_filename_component(_proc_dir "${MAIN_FILE}" DIRECTORY)
      if(EXISTS "${_proc_dir}/movenode.cc")
        list(APPEND SRCS "${_proc_dir}/movenode.cc")
      endif()
    endif()
    if(TARGET_NAME STREQUAL "plot_ugrid")
      # plot_ugrid main file name stem differs from implementation file providing InitThisProcess
      # Include the legacy implementation source so the symbol resolves.
      get_filename_component(_proc_dir "${MAIN_FILE}" DIRECTORY)
      if(EXISTS "${_proc_dir}/plot_misorient_unodebound.elle.cc")
        list(APPEND SRCS "${_proc_dir}/plot_misorient_unodebound.elle.cc")
      endif()
    endif()
    add_executable(${TARGET_NAME} ${SRCS})
    # Legacy Fortran source compatibility: many old fixed-form constructs and
    # deleted features (arithmetic IF, shared DO labels) appear in some
    # processes (e.g. tbh). Modern compilers warn but should still compile.
    # Provide permissive flags only for targets that actually have .f sources.
    foreach(_s IN LISTS SRCS)
      if(_s MATCHES "\\.f$")
        # Use target-specific compile options instead of global flags so newer
        # C++ sources remain on strict warning settings.
        # -std=legacy : accept deleted/obsolescent features (gfortran alias)
        # -ffixed-line-length-none : avoid truncation of legacy 72/80 cols
        # -fallow-argument-mismatch : relax old interface checks
        # (Some flags may be ignored by non-GNU compilers; that's acceptable.)
        target_compile_options(${TARGET_NAME} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy -ffixed-line-length-none -fallow-argument-mismatch>)
        break()
      endif()
    endforeach()
    if(TARGET_NAME STREQUAL "gbm")
      target_compile_definitions(${TARGET_NAME} PRIVATE XY=1)
    endif()
    target_include_directories(${TARGET_NAME} PRIVATE
      "${ELLE_BASECODE_DIR}"
      "${ELLE_MATRIX_DIR}"
      "${ELLE_STATSCODE_DIR}"
      "${ELLE_GPC_DIR}"
      "${CMAKE_SOURCE_DIR}/include")
    # Ensure stub headers (iphreeqc) are reachable for legacy processes using raw include
    target_include_directories(${TARGET_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/include")
    target_link_libraries(${TARGET_NAME} PRIVATE elle_core)
    if(UNIX)
      target_link_libraries(${TARGET_NAME} PRIVATE m)
    endif()
    # Link OpenMP if available (phreeqc.cc and similar use omp_get_max_threads)
    if(TARGET OpenMP::OpenMP_CXX)
      target_link_libraries(${TARGET_NAME} PRIVATE OpenMP::OpenMP_CXX)
      target_compile_definitions(${TARGET_NAME} PRIVATE ELLE_USE_OPENMP)
    endif()
    set_target_properties(${TARGET_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/processes"
      OUTPUT_NAME ${TARGET_NAME})
    list(APPEND ADDED_PROCESS_TARGETS ${TARGET_NAME})
  endforeach()
  if(NEEDS_FORTRAN)
    enable_language(Fortran)
  endif()
  if(ADDED_PROCESS_TARGETS)
    message(STATUS "Selected process executables: ${ADDED_PROCESS_TARGETS}")
    install(TARGETS ${ADDED_PROCESS_TARGETS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/processes)
  else()
    if(ELLE_PROCESS AND NOT ELLE_BUILD_ALL_PROCESSES)
      message(WARNING "Requested process '${ELLE_PROCESS}' not found (no matching *.main.cc)")
    endif()
  endif()
endif()

install(TARGETS latte_clean_old RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# ------------------------------------------------------------
# Process: latte (newer process) — headless and optional GUI
# ------------------------------------------------------------
set(LATTE_DIR "${CMAKE_SOURCE_DIR}/elle/elle/processes/latte")

set(LATTE_SOURCES
  "${LATTE_DIR}/elle.main.cc"
  "${LATTE_DIR}/experiment.cc"
  "${LATTE_DIR}/lattice.cc"
  "${LATTE_DIR}/particle.cc"
  "${LATTE_DIR}/phase_base.cc"
  "${LATTE_DIR}/phase_lattice.cc"
  "${LATTE_DIR}/fluid_lattice.cc"
  "${LATTE_DIR}/heat_lattice.cc"
  "${LATTE_DIR}/pressure_lattice.cc"
  "${LATTE_DIR}/graingrowth.cc"
  "${LATTE_DIR}/min_trans.cc"
  "${LATTE_DIR}/phreeqc.cc"
)

add_executable(latte ${LATTE_SOURCES})
target_include_directories(latte PRIVATE
  "${LATTE_DIR}"
  "${ELLE_BASECODE_DIR}"
  "${ELLE_MATRIX_DIR}"
  "${ELLE_STATSCODE_DIR}"
  "${ELLE_GPC_DIR}"
  "${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/include" # ensure include root is present
)
# Suppress deprecation error from boost::progress/timer (legacy code depends on it)
target_compile_definitions(latte PRIVATE LATTE_DATA_DIR="${CMAKE_SOURCE_DIR}/elle/elle/processes/latte")
target_link_libraries(latte PRIVATE elle_core Boost::boost)
if(TARGET OpenMP::OpenMP_CXX)
  target_link_libraries(latte PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(latte PRIVATE ELLE_USE_OPENMP)
endif()
if(IPhreeqc_FOUND)
  target_link_libraries(latte PRIVATE IPhreeqc::IPhreeqc)
endif()
if(UNIX)
  target_link_libraries(latte PRIVATE m)
endif()
install(TARGETS latte RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# GUI variant when available
if(ELLE_BUILD_GUI AND wxWidgets_FOUND AND ELLE_GUI_SOURCES)
  set(LATTE_GUI_SOURCES ${LATTE_SOURCES})
  # Exclude the CLI main when building the GUI app; wxstart.cc provides main
  list(FILTER LATTE_GUI_SOURCES EXCLUDE REGEX ".*/elle\\.main\\.cc$")
  # nodisp stubs are removed from build; GUI variant uses native display code
  add_executable(latte_gui ${LATTE_GUI_SOURCES} ${ELLE_GUI_SOURCES})
  target_compile_definitions(latte_gui PRIVATE BOOST_TIMER_ENABLE_DEPRECATED=1)
  # Provide data dir for runtime resource lookup as well for GUI
  target_compile_definitions(latte_gui PRIVATE LATTE_DATA_DIR="${CMAKE_SOURCE_DIR}/elle/elle/processes/latte")
  target_include_directories(latte_gui PRIVATE
    "${LATTE_DIR}"
    "${ELLE_BASECODE_DIR}"
    "${ELLE_MATRIX_DIR}"
    "${ELLE_STATSCODE_DIR}"
    "${ELLE_GPC_DIR}"
  "${CMAKE_SOURCE_DIR}/include")
  target_link_libraries(latte_gui PRIVATE elle_core_gui wxWidgets::wxWidgets)
  if(IPhreeqc_FOUND)
    target_link_libraries(latte_gui PRIVATE IPhreeqc::IPhreeqc)
  endif()
  if(UNIX)
    target_link_libraries(latte_gui PRIVATE m)
  endif()
  install(TARGETS latte_gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
