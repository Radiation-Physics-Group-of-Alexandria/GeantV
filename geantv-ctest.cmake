####################################################################
# Before run should be exported next variables:
# $CTEST_BUILD_OPTIONS // CMake flags for Geant-V build
# $CMAKE_SOURCE_DIR    // CMake source directory
# $CMAKE_BINARY_DIR    // CMake binary directory
# $CMAKE_BUILD_TYPE    // CMake build type: Debug, Release
# $CMAKE_INSTALL_PREFIX // Installation prefix for CMake (Jenkins trigger)
# CC and CXX (In Jenkins this step has been done authomaticly)
# export $LD_LIBRARY_PATH=$WORKSPACE/lib:$LD_LIBRARY_PATH (for GeantV libraries)
# Enviroment for name of build for CERN CDash:
# $LABEL                // Name of node (Jenkins trigger)
# Name of $BACKEND     // Backend for Geant-V (VecGeom/ROOT or CUDA)

cmake_minimum_required(VERSION 2.8)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS "1000")
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "1000")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "50000")
###################################################################
macro(CheckExitCode)
  if(NOT ${ExitCode} EQUAL 0)
    return(${ExitCode})
 endif()
endmacro(CheckExitCode)

####################################################################
# Build name settings
find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)

if(DEFINED ENV{LABEL})
  if (DEFINED ENV{BACKEND})
  set(CTEST_BUILD_NAME "${osname}-${cpu}-$ENV{LABEL}-$ENV{BACKEND}-$ENV{CMAKE_BUILD_TYPE}")
  endif()
else()
  set(CTEST_BUILD_NAME "${osname}-${cpu}-$ENV{LABEL}-$ENV{CMAKE_BUILD_TYPE}")
endif()
message("CTEST name: ${CTEST_BUILD_NAME}")

find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)
IF(NOT DEFINED CTEST_SITE)
  SET(CTEST_SITE "${HOSTNAME}")
ENDIF(NOT DEFINED CTEST_SITE)

#######################################################
# Build dashboard model setup

#------------ Cdash---------------------------------
if("$ENV{MODE}" STREQUAL "")
  set(CTEST_MODE Experimental)
else()
  set(CTEST_MODE "$ENV{MODE}")
endif()

if(${CTEST_MODE} MATCHES nightly)
  SET(MODEL Nightly)
elseif(${CTEST_MODE} MATCHES continuous)
  SET(MODEL Continuous)
elseif(${CTEST_MODE} MATCHES memory)
  SET(MODEL NightlyMemoryCheck)
else()
  SET(MODEL Experimental)
endif()

find_program(CTEST_COMMAND_BIN NAMES ctest)
SET (CTEST_COMMAND
    "$CTEST_COMMAND_BIN -D ${MODEL}")

#######################################################
if(MODEL STREQUAL NightlyMemoryCheck)
  set(WITH_MEMCHECK TRUE)
#  find_package(Valgrind)
#  find_package(ROOT REQUIRED)
  if(ROOT)
    set(MEMORYCHECK_SUPPRESSIONS_FILE "${ROOTSYS}/etc/valgrind-root.supp")
  endif()
  set(CTEST_MEMORYCHECK_COMMAND ${VALGRIND_PROGRAM})
else()
  set(WITH_MEMCHECK FALSE)
endif()

######################################################
set(WITH_COVERAGE FALSE)

#######################################################
# CTest/CMake settings

set(CTEST_TEST_TIMEOUT 3600)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS "2000")
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "2000")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "5000")
set(CTEST_BUILD_CONFIGURATION "$ENV{CMAKE_BUILD_TYPE}")
set(CMAKE_INSTALL_PREFIX "$ENV{CMAKE_INSTALL_PREFIX}")
set(CTEST_SOURCE_DIRECTORY "$ENV{CMAKE_SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY "$ENV{CMAKE_BINARY_DIR}")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_OPTIONS "$ENV{CTEST_BUILD_OPTIONS}")
set(CTEST_CONFIGURE_COMMAND "\"${CMAKE_COMMAND}\" \"-G${CTEST_CMAKE_GENERATOR}\" -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE=${CTEST_BUILD_CONFIGURATION} ${CTEST_BUILD_OPTIONS} ${CTEST_SOURCE_DIRECTORY}")
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

#########################################################
# git command configuration

find_program(CTEST_GIT_COMMAND NAMES git)
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone https://gitlab.cern.ch/GeantV/geant.git ${CTEST_SOURCE_DIRECTORY}")
endif()
set(CTEST_GIT_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
if(${MODEL} MATCHES Nightly OR Experimental)
    if(NOT "$ENV{GIT_COMMIT}" STREQUAL "")
        set(CTEST_CHECKOUT_COMMAND "cmake -E chdir ${CTEST_SOURCE_DIRECTORY} ${CTEST_GIT_COMMAND} checkout -f $ENV{GIT_PREVIOUS_COMMIT}")
        set(CTEST_GIT_UPDATE_CUSTOM  ${CTEST_GIT_COMMAND} checkout -f $ENV{GIT_COMMIT})
    endif()
else()
    if(NOT "$ENV{GIT_COMMIT}" STREQUAL "")
        set(CTEST_CHECKOUT_COMMAND "cmake -E chdir ${CTEST_SOURCE_DIRECTORY} ${CTEST_GIT_COMMAND} checkout -f $ENV{GIT_COMMIT}")
        set(CTEST_GIT_UPDATE_CUSTOM  ${CTEST_GIT_COMMAND} checkout -f $ENV{GIT_COMMIT})
    endif()
endif()

#########################################################
## Output language
set($ENV{LC_MESSAGES}  "en_EN")

#########################################################
# Use multiple CPU cores to build
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
  if(NOT WIN32)
#    set(CTEST_BUILD_FLAGS -j${N})

  endif(NOT WIN32)
  set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()

##########################################################
# Print summary information.
foreach(v
    CTEST_SITE
    CTEST_BUILD_NAME
    CTEST_SOURCE_DIRECTORY
    CTEST_BINARY_DIRECTORY
    CTEST_CMAKE_GENERATOR
    CTEST_BUILD_CONFIGURATION
    CTEST_GIT_COMMAND
    CTEST_CONFIGURE_COMMAND
    CTEST_SCRIPT_DIRECTORY
    CTEST_BUILD_FLAGS
    WITH_MEMCHECK
    WITH_COVERAGE
  )
  set(vars "${vars}  ${v}=[${${v}}]\n")
endforeach(v)
message("Dashboard script configuration (check if everything is declared correctly):\n${vars}\n")

#######################################################
# Test custom update with a dashboard script.
message("Running CTest Dashboard Script (custom update)...")
include("${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake")

#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

ctest_start(${MODEL} TRACK ${MODEL})

ctest_update(SOURCE ${CTEST_SOURCE_DIRECTORY})
message("Updated.")

ctest_configure(SOURCE "${CTEST_SOURCE_DIRECTORY}" BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
message("Configured.")
ctest_submit(PARTS Update Configure Notes)

ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
message("Built.")
ctest_submit(PARTS Build)

message(" -- Install ${MODEL} - ${CTEST_BUILD_NAME} --")
execute_process(COMMAND make install  WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}  RESULT_VARIABLE ExitCode)
CheckExitCode()

ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
message("Tested.")
ctest_submit(PARTS Test)

if(${MODEL} MATCHES NightlyMemoryCheck)
  ctest_submit(PARTS MemCheck)
endif()

message("DONE:CTestScript")
