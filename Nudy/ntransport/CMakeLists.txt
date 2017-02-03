cmake_minimum_required(VERSION 2.8)
project(NudyTest)

file(GLOB sources src/*.cxx)
file(GLOB headers inc/.hh)

include_directories(inc)

include(${Geant4_USE_FILE})

include_directories(include 
  ${Geant4_INCLUDE_DIR})

#add_executable(nfission nfission.cpp)
#target_link_libraries(nfission Nudy)


add_executable(ninelastic ninelastic.cpp)
target_link_libraries(ninelastic Nudy)

#add_executable(ncapture ncapture.cpp)
#target_link_libraries(ncapture Nudy)


#add_executable(nelastic nelastic.cpp)
#target_link_libraries(nelastic Nudy)


#add_executable(nthermal nthermal.cpp)
#target_link_libraries(nthermal Nudy)
