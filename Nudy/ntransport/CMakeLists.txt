cmake_minimum_required(VERSION 2.8)
project(NudyTest)

file(GLOB sources src/*.cxx)
file(GLOB headers inc/.hh)

include_directories(inc)

include(${Geant4_USE_FILE})

include_directories(include 
  ${Geant4_INCLUDE_DIR})
#add_executable(dumy2 dumy2.cpp)
#target_link_libraries(dumy2 Nudy)


#add_executable(fission fission.cpp)
#target_link_libraries(fission Nudy)


add_executable(inelastic inelastic.cpp)
target_link_libraries(inelastic Nudy)


#add_executable(thermal thermal.cpp)
#target_link_libraries(thermal Nudy)




#add_executable(thermal thermal.cpp)
#target_link_libraries(thermal Nudy)
