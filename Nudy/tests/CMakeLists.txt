file (GLOB sources src/*.cxx)
file(GLOB headers inc/.h)

include_directories(inc)

include(${Geant4_USE_FILE})

include_directories(include ${Geant4_INCLUDE_DIR})

add_executable(nelastic src/nelastic.cxx)

#add_library(Nudy SHARED ${sources} NudyDict.cxx)
target_link_libraries(nelastic -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Nudy -L${ROOTSYS}/lib MathMore RIO Geom Graf Gpad Physics EG ${VC_LIBRARIES})

target_link_libraries(nelastic Nudy)
