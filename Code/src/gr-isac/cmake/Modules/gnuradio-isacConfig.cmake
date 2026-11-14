find_package(PkgConfig)

PKG_CHECK_MODULES(PC_GR_ISAC gnuradio-isac)

FIND_PATH(
    GR_ISAC_INCLUDE_DIRS
    NAMES gnuradio/isac/api.h
    HINTS $ENV{ISAC_DIR}/include
        ${PC_ISAC_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    GR_ISAC_LIBRARIES
    NAMES gnuradio-isac
    HINTS $ENV{ISAC_DIR}/lib
        ${PC_ISAC_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

include("${CMAKE_CURRENT_LIST_DIR}/gnuradio-isacTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GR_ISAC DEFAULT_MSG GR_ISAC_LIBRARIES GR_ISAC_INCLUDE_DIRS)
MARK_AS_ADVANCED(GR_ISAC_LIBRARIES GR_ISAC_INCLUDE_DIRS)
