# Copyright 2020 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

########################################################################
# Check if there is C++ code at all
########################################################################
if(NOT isac_sources)
    message(STATUS "No C++ sources... skipping python bindings")
    return()
endif(NOT isac_sources)

########################################################################
# Check for pygccxml
########################################################################
gr_python_check_module_raw("pygccxml" "import pygccxml" PYGCCXML_FOUND)

include(GrPybind)

########################################################################
# Python Bindings
########################################################################

list(APPEND isac_python_files python_bindings.cc)

gr_pybind_make_oot(isac ../../.. gr::isac "${isac_python_files}")

# copy bindings extension for use in QA test module
add_custom_command(
    TARGET isac_python
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:isac_python>
            ${PROJECT_BINARY_DIR}/test_modules/gnuradio/isac/)

install(
    TARGETS isac_python
    DESTINATION ${GR_PYTHON_DIR}/gnuradio/isac
    COMPONENT pythonapi)
