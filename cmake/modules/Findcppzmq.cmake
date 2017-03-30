# Copied from: https://github.com/zeromq/lzmq/blob/master/cmake/FindZMQ.cmake
#
# - Try to find ZMQ
# Once done this will define
# ZMQ_FOUND - System has ZMQ
# ZMQ_INCLUDE_DIRS - The ZMQ include directories
# ZMQ_LIBRARIES - The libraries needed to use ZMQ
# ZMQ_DEFINITIONS - Compiler switches required for using ZMQ

find_path(cppzmq_INCLUDE_DIR zmq.h)
find_library(cppzmq_LIBRARY NAMES zmq)

set(cppzmq_LIBRARIES ${ZMQ_LIBRARY})
set(cppzmq_INCLUDE_DIRS ${ZMQ_INCLUDE_DIR})

include (FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set ZMQ_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(cppzmq DEFAULT_MSG cppzmq_LIBRARY cppzmq_INCLUDE_DIR )

