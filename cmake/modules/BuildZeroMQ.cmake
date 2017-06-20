set(ZMQ_CMAKE_FILE_CONTENT
        "
cmake_minimum_required(VERSION 3.4)
include(ExternalProject)
ExternalProject_Add(Zeromq
  GIT_REPOSITORY https://github.com/zeromq/libzmq.git

  LOG_DOWNLOAD 1 LOG_CONFIGURE 1 LOG_BUILD 1 LOG_INSTALL 1
  CMAKE_ARGS  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
               -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
               -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install)

ExternalProject_Add(cppzmq
  DEPENDS Zeromq
  GIT_REPOSITORY https://github.com/zeromq/cppzmq.git

  LOG_DOWNLOAD 1 LOG_CONFIGURE 1 LOG_BUILD 1 LOG_INSTALL 1
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
               -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
               -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install)"
        )

FILE(WRITE "${CMAKE_BINARY_DIR}/external/CMakeLists.txt" "${ZMQ_CMAKE_FILE_CONTENT}")

execute_process(
        COMMAND ${CMAKE_COMMAND} ${CMAKE_BINARY_DIR}/external/
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/
)
execute_process(
        COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/
)

unset(ZMQ_CMAKE_FILE_CONTENT)
